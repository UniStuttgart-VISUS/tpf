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
        /// <param name="data_time_range">Available time range</param>
        /// <param name="request">Original request that will be modified for requesting different time frames</param>
        /// <param name="grid_alg">Algorithm producing the velocity field and additional information</param>
        /// <param name="droplets_alg">Algorithm producing droplet data</param>
        /// <param name="get_array_name">Function to query the respective array names</param>
        /// <param name="array_selection">Selection of property arrays to interpolate and store at particle positions</param>
        /// <param name="initial_time">Initial time for integration</param>
        /// <param name="timestep">Timestep size (delta t)</param>
        data_handler(const std::array<double, 2>& data_time_range, vtkInformation* request,
            vtkAlgorithm* grid_alg, vtkAlgorithm* droplets_alg,
            std::function<std::string(int, int)> get_array_name, vtkDataArraySelection* array_selection,
            const double initial_time, const double timestep)
        {
            this->time_offset = this->original_time = initial_time;
            this->data_time_range = data_time_range;

            this->request = request;
            this->grid_alg = grid_alg;
            this->droplets_alg = droplets_alg;

            this->get_array_name = get_array_name;

            this->array_selection = array_selection;

            this->timestep = timestep;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, sample time, velocities, translation, angular velocity, droplet barycenter,
        ///  initial translation, initial angular velocity, validity,
        ///  fields to interpolate and store at the particle positions]</returns>
        virtual std::tuple<
            float_t, float_t,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::vector<std::tuple<std::string, std::size_t, tpf::policies::interpolatable_base<Eigen::Matrix<float_t, 3, 1>>*>>> operator()() override
        {
            // Get data
            if ((this->time_offset >= this->data_time_range[0] && this->time_offset <= this->data_time_range[1]))
            {
                // Load requested time step
                this->request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->time_offset);

                this->grid_alg->Update(this->request);

                if (this->droplets_alg != nullptr)
                {
                    this->droplets_alg->Update(this->request);
                }

                // Load input grids and droplets
                auto in_grid = vtkRectilinearGrid::SafeDownCast(this->grid_alg->GetOutputDataObject(0));

                this->velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(2, 0));

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
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_initial_translation;
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_initial_angular_velocity;
                std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)> is_valid;

                if (this->droplets_alg != nullptr && vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0)) != nullptr)
                {
                    auto in_droplets = vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0));

                    this->droplet_grid = std::make_shared<tpf::data::grid<long long, float_t, 3, 1>>(
                        tpf::vtk::get_grid<long long, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(1, 0)));
                    const auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
                        tpf::data::data_information<float_t, 3>{ get_array_name(4, 1), tpf::data::topology_t::POINT_DATA },
                        tpf::data::data_information<float_t, 3>{ get_array_name(5, 1), tpf::data::topology_t::POINT_DATA });

                    // Store angular velocities and droplet velocities
                    this->droplets = droplets.get_geometry();
                    this->droplet_velocity = droplets.template get_point_data_as<float_t, 3>(get_array_name(4, 1));
                    this->angular_velocity = droplets.template get_point_data_as<float_t, 3>(get_array_name(5, 1));

                    if (this->original_time == this->time_offset)
                    {
                        this->previous_droplet_grid = this->droplet_grid;
                        this->initial_droplet_velocity = this->droplet_velocity;
                        this->initial_angular_velocity = this->angular_velocity;

                        for (std::size_t i = 0; i < this->droplets.size(); ++i)
                        {
                            this->initial_droplet[i] = i;
                        }
                    }
                    else
                    {
                        // Advect all droplet barycenters backwards and look up previous droplet ID
                        std::map<long long, long long> new_mapping;

                        for (std::size_t i = 0; i < this->droplets.size(); ++i)
                        {
                            const auto barycenter = this->droplets[i]->get_points().front();
                            const auto advected = barycenter - this->timestep * this->droplet_velocity->at(i);

                            const auto previous_cell = this->previous_droplet_grid->find_cell(advected);

                            if (previous_cell)
                            {
                                const auto previous_id = (*this->previous_droplet_grid)(*previous_cell);

                                if (this->initial_droplet.find(previous_id) != this->initial_droplet.end())
                                {
                                    new_mapping[i] = this->initial_droplet.at(previous_id);
                                }
                            }
                        }

                        std::swap(this->initial_droplet, new_mapping);

                        this->previous_droplet_grid = this->droplet_grid;
                    }

                    // Create look-up functions
                    get_translation = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            const auto id = (*this->droplet_grid)(*cell);

                            return this->droplet_velocity->at(id);
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_angular_velocity = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            const auto id = (*this->droplet_grid)(*cell);

                            return this->angular_velocity->at(id);
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_barycenter = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            const auto id = (*this->droplet_grid)(*cell);

                            return this->droplets[id]->get_points().front();
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_initial_translation = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            const auto id = (*this->droplet_grid)(*cell);

                            return this->initial_droplet_velocity->at(this->initial_droplet.at(id));
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_initial_angular_velocity = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            const auto id = (*this->droplet_grid)(*cell);

                            return this->initial_angular_velocity->at(this->initial_droplet.at(id));
                        }

                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    is_valid = [this](const Eigen::Matrix<float_t, 3, 1>& position)
                    {
                        const auto cell = this->droplet_grid->find_cell(position);

                        if (cell)
                        {
                            return (*this->droplet_grid)(*cell) >= 0;
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

                    get_initial_translation = [](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    get_initial_angular_velocity = [](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                    {
                        return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                    };

                    is_valid = [this](const Eigen::Matrix<float_t, 3, 1>& position) -> bool
                    {
                        return static_cast<bool>(this->velocity_grid.find_cell(position));
                    };
                }

                // Get time information and set new timestep for next iteration
                if (this->data_time_range[0] != this->data_time_range[1])
                {
                    this->time_offset += this->timestep;
                }

                const auto sample_time = this->grid_alg->GetOutputDataObject(0)->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());

                // Return [Time step delta, velocities, translation, angular velocity, barycenter, initial translation,
                //  initial angular velocity, validity, fields to interpolate and store at the particle positions]
                return std::make_tuple(this->timestep, sample_time,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*>(&this->velocity_grid),
                    get_translation, get_angular_velocity, get_barycenter, get_initial_translation, get_initial_angular_velocity, is_valid, property_grids);
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

            request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->original_time);

            this->grid_alg->Update(request);

            if (this->droplets_alg != nullptr)
            {
                this->droplets_alg->Update(request);
            }
        }

    private:
        /// Original timestep
        double original_time;

        /// Last requested timestep
        double time_offset;

        /// Available time range
        std::array<double, 2> data_time_range;

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

        /// Timestep information
        double timestep;

        /// Data of current timestep
        tpf::data::grid<float_t, float_t, 3, 3> velocity_grid;
        std::shared_ptr<tpf::data::grid<long long, float_t, 3, 1>> droplet_grid, previous_droplet_grid;

        std::vector<std::shared_ptr<tpf::policies::interpolatable_base<Eigen::Matrix<float_t, 3, 1>>>> property_grids;

        /// Droplets, their angular velocity and velocity
        std::vector<std::shared_ptr<tpf::geometry::geometric_object<float_t>>> droplets;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> angular_velocity;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> droplet_velocity;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> initial_angular_velocity;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> initial_droplet_velocity;

        /// Map to initial droplet index
        std::map<long long, long long> initial_droplet;
    };
}

vtkStandardNewMacro(tpf_flow_field);

tpf_flow_field::tpf_flow_field() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(3);
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
    else if (port == 2)
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
        auto in_droplets = vtkPolyData::GetData(input_vector[2]);

        const auto get_array_name = [this, &in_grid, &in_droplets](const int index, const int data_index = 0) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, data_index == 0 ? static_cast<vtkDataSet*>(in_grid) : static_cast<vtkDataSet*>(in_droplets));

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")"));
        };

        std::array<double, 2> data_time_range;
        input_vector[0]->GetInformationObject(0)->Get(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), data_time_range.data());
        const auto current_timestep = in_grid->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());

        const std::array<double, 2> integration_range = {
            std::min(this->TimeRange[0], this->TimeRange[1]),
            std::max(this->TimeRange[0], this->TimeRange[1]) };

        const std::size_t num_advections =
            (static_cast<tpf::modules::flow_field_aux::method_t>(this->Method) == tpf::modules::flow_field_aux::method_t::STREAM)
            ? static_cast<std::size_t>(this->NumAdvections)
            : static_cast<std::size_t>(std::ceil((integration_range[1] - integration_range[0]) / this->FixedTimeStep));

        data_handler<float_t> call_back_loader(data_time_range, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(2, 0),
            get_array_name, this->array_selection, ((this->FixedTimeStep > 0.0) ? integration_range[0] : integration_range[1]), this->FixedTimeStep);

        // Create seed
        tpf::data::polydata<float_t> seed;

        if (this->SeedInCells != 0)
        {
            auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, get_array_name(0, 0));

            for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                    {
                        const tpf::data::coords3_t coords(x, y, z);

                        if ((this->SeedCellType == 0 && vof(coords) > 0.0) ||
                            (this->SeedCellType == 1 && vof(coords) > 0.0 && vof(coords) < 1.0) ||
                            (this->SeedCellType == 2 && vof(coords) == 1.0))
                        {
                            seed.insert(std::make_shared<tpf::geometry::point<float_t>>(vof.get_cell_coordinates(coords)));
                        }
                    }
                }
            }
        }
        else
        {
            auto in_seed = vtkPolyData::GetData(input_vector[1]);

            if (in_seed != nullptr)
            {
                seed = tpf::vtk::get_polydata<float_t>(in_seed);
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("No input seed available."));
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
            num_advections, static_cast<tpf::modules::flow_field_aux::time_dependency_t>(this->TimeDependency),
            this->KeepTranslation == 1, this->KeepRotation == 1);

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
