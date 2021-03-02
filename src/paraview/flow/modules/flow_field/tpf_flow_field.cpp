#include "tpf_flow_field.h"

#include "tpf_modules/flow_field/tpf_module_flow_field.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/policies/tpf_interpolatable.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkAlgorithm.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkFieldData.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
    /// Call back class for next time step
    /// <template name="float_t">Floating point type</template>
    template <typename float_t>
    class data_handler : public tpf::modules::flow_field_aux::request_frame_call_back<float_t, Eigen::Matrix<float_t, 3, 1>>
    {
    public:
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="current_timestep">Current time step this plugin was executed in</param>
        /// <param name="timesteps">Available time steps</param>
        /// <param name="request">Original request that will be modified for requesting different time frames</param>
        /// <param name="grid_alg">Algorithm producing the velocity field and additional information</param>
        /// <param name="droplets_alg">Algorithm producing droplet data</param>
        /// <param name="get_array_name">Function to query the respective array names</param>
        /// <param name="fixed_timestep">Optional fixed timestep</param>
        /// <param name="time_from_data">Extract timestep information from data</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps,
            vtkInformation* request, vtkAlgorithm* grid_alg, vtkAlgorithm* droplets_alg, std::function<std::string(int, int)> get_array_name,
            vtkDataArraySelection* array_selection, std::optional<float_t> fixed_timestep = std::nullopt, const bool time_from_data = false)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->grid_alg = grid_alg;
            this->droplets_alg = droplets_alg;

            this->get_array_name = get_array_name;

            this->array_selection = array_selection;

            this->fixed_timestep = fixed_timestep;
            this->time_from_data = time_from_data;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, velocities, global velocity parts, translation,
        ///  angular velocity, validity, fields to interpolate and store at the particle positions]</returns>
        virtual std::tuple<
            float_t,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::vector<std::tuple<std::string, std::size_t, tpf::policies::interpolatable_base<Eigen::Matrix<float_t, 3, 1>>*>>> operator()() override
        {
            // Get data
            if (this->timesteps.size() > this->time_offset || this->original_time == this->time_offset)
            {
                // Load input grids and droplets
                auto in_grid = vtkRectilinearGrid::SafeDownCast(this->grid_alg->GetOutputDataObject(0));

                this->droplet_grid = tpf::vtk::get_grid<long long, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(1, 0));
                this->velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(2, 0));
                this->global_velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(3, 0));

                // Get property arrays
                std::vector<vtkDataArray*> property_arrays;

                for (int index = 0; index < in_grid->GetPointData()->GetNumberOfArrays(); ++index)
                {
                    auto in_array = in_grid->GetPointData()->GetArray(index);

                    if (in_array && in_array->GetName())
                    {
                        if (this->array_selection->ArrayIsEnabled(in_array->GetName()))
                        {
                            property_arrays.push_back(in_array);
                        }
                    }
                }

                // Create grids for property arrays
                this->property_grids.resize(property_arrays.size());

                std::vector<std::tuple<std::string, std::size_t, tpf::policies::interpolatable_base<Eigen::Matrix<float_t, 3, 1>>*>>
                    property_grids(property_arrays.size());

                for (std::size_t i = 0; i < property_arrays.size(); ++i)
                {
                    if (property_arrays[i]->GetNumberOfComponents() == 3)
                    {
                        this->property_grids[i] = std::make_shared<tpf::data::grid<double, float_t, 3, 3>>(
                            tpf::vtk::get_grid<double, float_t, 3, 3>(in_grid, tpf::data::topology_t::POINT_DATA, property_arrays[i]->GetName()));

                        property_grids[i] = std::make_tuple(std::string(property_arrays[i]->GetName()), 3, this->property_grids[i].get());
                    }
                    else if (property_arrays[i]->GetNumberOfComponents() == 1)
                    {
                        this->property_grids[i] = std::make_shared<tpf::data::grid<double, float_t, 3, 1>>(
                            tpf::vtk::get_grid<double, float_t, 3, 1>(in_grid, tpf::data::topology_t::POINT_DATA, property_arrays[i]->GetName()));

                        property_grids[i] = std::make_tuple(std::string(property_arrays[i]->GetName()), 1, this->property_grids[i].get());
                    }
                    else
                    {
                        tpf::log::info_message(__tpf_warning_message("Only property arrays with 1 or 3 components are supported."));
                    }
                }

                // Create look-up functions
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_translation;
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_angular_velocity;
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_barycenter;
                std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)> is_valid;

                if (this->droplets_alg != nullptr && vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0)) != nullptr)
                {
                    auto in_droplets = vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0));

                    const auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
                        tpf::data::data_information<float_t, 3>{ get_array_name(4, 1), tpf::data::topology_t::POINT_DATA },
                        tpf::data::data_information<float_t, 3>{ get_array_name(5, 1), tpf::data::topology_t::POINT_DATA });

                    // Store angular velocities and droplet velocities
                    this->droplets = droplets.get_geometry();
                    this->droplet_velocity = droplets.template get_point_data_as<float_t, 3>(get_array_name(4, 1));
                    this->angular_velocity = droplets.template get_point_data_as<float_t, 3>(get_array_name(5, 1));

                    // Create look-up functions
                    get_translation = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid.find_cell(position);

                        if (cell)
                        {
                            const auto id = this->droplet_grid(*cell);

                            return this->droplet_velocity->at(id);
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_angular_velocity = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid.find_cell(position);

                        if (cell)
                        {
                            const auto id = this->droplet_grid(*cell);

                            return this->angular_velocity->at(id);
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_barycenter = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid.find_cell(position);

                        if (cell)
                        {
                            const auto id = this->droplet_grid(*cell);

                            return this->droplets[id]->get_points().front();
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    is_valid = [this](const Eigen::Matrix<float_t, 3, 1>& position)
                    {
                        const auto cell = this->droplet_grid.find_cell(position);

                        if (cell)
                        {
                            return this->droplet_grid(*cell) >= 0;
                        }

                        return false;
                    };
                }
                else
                {
                    get_translation = [](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_angular_velocity = [](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_barycenter = [](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    is_valid = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> bool
                    {
                        return static_cast<bool>(this->droplet_grid.find_cell(position));
                    };
                }

                // Compute time step
                float_t timestep_delta;

                if (this->fixed_timestep)
                {
                    timestep_delta = *this->fixed_timestep;
                }
                else if (this->time_from_data)
                {
                    const auto time_array = in_grid->GetFieldData()->GetArray(get_array_name(6, 0).c_str());

                    if (time_array == nullptr)
                    {
                        tpf::log::warning_message(__tpf_warning_message("No valid time step data found. Using 1.0 instead."));

                        timestep_delta = 1.0;
                    }
                    else
                    {
                        timestep_delta = time_array->GetComponent(0, 0);
                    }
                }
                else if (this->timesteps.empty() || this->timesteps.size() == 1)
                {
                    tpf::log::warning_message(__tpf_warning_message("No valid time step found. Using 1.0 instead."));

                    timestep_delta = 1.0;
                }
                else if (this->time_offset == 0)
                {
                    timestep_delta = this->timesteps[1] - this->timesteps[0];
                }
                else
                {
                    timestep_delta = this->timesteps[this->time_offset] - this->timesteps[this->time_offset - 1];
                }

                if (timestep_delta == 0.0)
                {
                    timestep_delta = static_cast<float_t>(1.0);
                }

                // Check if next time step is available
                ++this->time_offset;

                if (this->timesteps.size() > this->time_offset)
                {
                    // Request update
                    this->request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                    this->grid_alg->Update(this->request);

                    if (this->droplets_alg != nullptr)
                    {
                        this->droplets_alg->Update(this->request);
                    }
                }

                // Return [Time step delta, velocities, global velocity parts, translation, angular velocity, validity]
                return std::make_tuple(timestep_delta,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*>(&this->velocity_grid),
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*>(&this->global_velocity_grid),
                    get_translation, get_angular_velocity, get_barycenter, is_valid, property_grids);
            }

            throw std::exception();
        }

        /// <summary>
        /// Reset input algorithms to their beginning state
        /// </summary>
        virtual void reset() override
        {
            // Request update to original time step
            this->time_offset = this->original_time;

            if (!this->timesteps.empty())
            {
                request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);
            }

            this->grid_alg->Update(request);

            if (this->droplets_alg != nullptr)
            {
                this->droplets_alg->Update(request);
            }
        }

    private:
        /// Original timestep
        std::size_t original_time;

        /// Last requested timestep
        std::size_t time_offset;

        /// Available timesteps
        std::vector<float_t> timesteps;

        /// Request
        vtkInformation* request;

        /// Grid algorithm
        vtkAlgorithm* grid_alg;

        /// Droplet algorithm
        vtkAlgorithm* droplets_alg;

        /// Function to retrieve array names
        std::function<std::string(int, int)> get_array_name;

        /// Property arrays
        vtkDataArraySelection* array_selection;

        /// Fixed timestep
        std::optional<float_t> fixed_timestep;
        bool time_from_data;

        /// Data of current timestep
        tpf::data::grid<float_t, float_t, 3, 3> velocity_grid, global_velocity_grid;
        tpf::data::grid<long long, float_t, 3, 1> droplet_grid;

        std::vector<std::shared_ptr<tpf::policies::interpolatable_base<Eigen::Matrix<float_t, 3, 1>>>> property_grids;

        /// Droplets, their angular velocity and velocity
        std::vector<std::shared_ptr<tpf::geometry::geometric_object<float_t>>> droplets;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> angular_velocity;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> droplet_velocity;
    };
}

vtkStandardNewMacro(tpf_flow_field);

tpf_flow_field::tpf_flow_field() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);

    this->array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);
}

tpf_flow_field::~tpf_flow_field() {}

vtkDataArraySelection* tpf_flow_field::GetArraySelection()
{
    return this->array_selection;
}

int tpf_flow_field::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }
    else if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
        return 1;
    }

    return 0;
}

int tpf_flow_field::RequestUpdateExtent(vtkInformation *vtkNotUsed(request), vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(), this->num_ghost_levels);

    return 1;
}

int tpf_flow_field::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);
        auto in_droplets = vtkPolyData::GetData(input_vector[1]);

        const auto get_array_name = [this, &in_grid, &in_droplets](const int index, const int data_index = 0) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, data_index == 0 ? static_cast<vtkDataSet*>(in_grid) : static_cast<vtkDataSet*>(in_droplets));

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")"));
        };

        auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(0, 0));

        const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_grid).second;
        const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

        data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(1, 0),
            get_array_name, this->array_selection, this->ForceFixedTimeStep ? std::make_optional(static_cast<float_t>(this->StreamTimeStep)) : std::nullopt,
            this->TimeStepFromData != 0);

        // Create seed
        tpf::data::polydata<float_t> seed;

        for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
        {
            for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
            {
                for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                {
                    const tpf::data::coords3_t coords(x, y, z);

                    if (vof(coords) > 0.0)
                    {
                        seed.insert(std::make_shared<tpf::geometry::point<float_t>>(vof.get_cell_coordinates(coords)));
                    }
                }
            }
        }

        // Create output data
        tpf::data::polydata<float_t> lines;

        // Set input, output and parameters
        using flow_t = tpf::modules::flow_field<float_t, Eigen::Matrix<float_t, 3, 1>>;

        flow_t flow_field_module;

        flow_field_module.set_input(seed);
        flow_field_module.set_output(lines);

        flow_field_module.set_parameters(static_cast<tpf::modules::flow_field_aux::method_t>(this->Method),
            static_cast<std::size_t>(this->NumAdvections));

        flow_field_module.set_callbacks(&call_back_loader);

        flow_field_module.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        tpf::vtk::set_polydata(output, lines,
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Advection)"), tpf::data::topology_t::CELL_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Distribution)"), tpf::data::topology_t::CELL_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Advection)"), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Distribution)"), tpf::data::topology_t::POINT_DATA });

        for (const auto& pd : lines.get_point_data())
        {
            if (pd->get_name().substr(0, 2).compare("ID"))
            {
                tpf::vtk::set_data<double>(output, tpf::data::topology_t::POINT_DATA, pd->get_name(),
                    pd->get_data_dynamic(), pd->get_num_components_dynamic());
            }
        }
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Flow Field' failed."));

        return 0;
    }

    return 1;
}
