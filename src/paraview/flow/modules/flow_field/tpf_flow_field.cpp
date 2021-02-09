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

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
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
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps,
            vtkInformation* request, vtkAlgorithm* grid_alg, vtkAlgorithm* droplets_alg)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->grid_alg = grid_alg;
            this->droplets_alg = droplets_alg;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, velocities, global velocity parts, translation, rotation, validity]</returns>
        virtual std::tuple<
            float_t,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)>,
            std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)>> operator()() override
        {
            // Get data
            if (this->timesteps.size() > this->time_offset || this->original_time == this->time_offset)
            {
                // Load input grids and droplets
                auto in_grid = vtkRectilinearGrid::SafeDownCast(this->grid_alg->GetOutputDataObject(0));
                auto in_droplets = vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0));

                this->droplet_grid = tpf::vtk::get_grid<long long, float_t, 3, 1>(in_grid,
                    tpf::data::topology_t::CELL_DATA, in_grid->GetPointData()->GetArray(get_array_name(1, 0).c_str()));
                this->velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid,
                    tpf::data::topology_t::CELL_DATA, in_grid->GetPointData()->GetArray(get_array_name(2, 0).c_str()));
                this->global_velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid,
                    tpf::data::topology_t::CELL_DATA, in_grid->GetPointData()->GetArray(get_array_name(3, 0).c_str()));

                auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
                    tpf::data::data_information<float_t, 3>{ get_array_name(4, 1), tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 3>{ get_array_name(5, 1), tpf::data::topology_t::POINT_DATA });

                // Store rotation axes and droplet velocities
                this->droplet_velocity = droplets.template get_point_data_as<float_t, 3>(get_array_name(4, 1));
                this->rotation_axis = droplets.template get_point_data_as<float_t, 3>(get_array_name(5, 1));

                // Create look-up functions
                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_translation =
                    [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                {
                    const auto cell = this->droplet_grid.find_cell(position);

                    if (cell)
                    {
                        const auto id = this->droplet_grid(*cell);

                        return this->droplet_velocity->at(id);
                    }

                    return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                };

                std::function<Eigen::Matrix<float_t, 3, 1>(const Eigen::Matrix<float_t, 3, 1>&)> get_rotation =
                    [this](const Eigen::Matrix<float_t, 3, 1>& position) -> Eigen::Matrix<float_t, 3, 1>
                {
                    const auto cell = this->droplet_grid.find_cell(position);

                    if (cell)
                    {
                        const auto id = this->droplet_grid(*cell);

                        return this->rotation_axis->at(id);
                    }

                    return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0);
                };

                std::function<bool(const Eigen::Matrix<float_t, 3, 1>&)> is_valid =
                    [this](const Eigen::Matrix<float_t, 3, 1>& position)
                {
                    const auto cell = this->droplet_grid.find_cell(position);

                    if (cell)
                    {
                        return this->droplet_grid(*cell) >= 0;
                    }

                    return false;
                };

                // Compute time step
                const float_t timestep_delta = this->fixed_timestep ? *fixed_timestep : (this->timesteps[this->time_offset] - this->timesteps[this->time_offset - 1]);

                // Check if next time step is available
                ++this->time_offset;

                if (this->timesteps.size() > this->time_offset)
                {
                    // Request update
                    this->request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                    this->grid_alg->Update(this->request);
                    this->droplets_alg->Update(this->request);
                }

                // Return [Time step delta, velocities, global velocity parts, translation, rotation, validity]
                return std::make_tuple(timestep_delta,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*>(&this->velocity_grid),
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>>*>(&this->global_velocity_grid),
                    get_translation, get_rotation, is_valid);
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

            request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

            this->grid_alg->Update(request);
            this->droplets_alg->Update(request);
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

        /// Fixed timestep
        std::optional<float_t> fixed_timestep;

        /// Data of current timestep
        tpf::data::grid<float_t, float_t, 3, 3> velocity_grid, global_velocity_grid;
        tpf::data::grid<long long, float_t, 3, 1> droplet_grid;

        /// Rotation and velocity
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> rotation_axis;
        std::shared_ptr<tpf::data::array<float_t, 3, 1>> droplet_velocity;
    };
}

vtkStandardNewMacro(tpf_flow_field);

tpf_flow_field::tpf_flow_field() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

tpf_flow_field::~tpf_flow_field() {}

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

        auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, in_grid->GetPointData()->GetArray(get_array_name(0, 0).c_str()));

        const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_grid).second;
        const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

        data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(1, 0));

        const auto initial_data = call_back_loader();
        const auto initial_velocities = std::get<1>(initial_data);
        const auto initial_global_velocities = std::get<2>(initial_data);
        const auto& initial_get_translation = std::get<3>(initial_data);
        const auto& initial_get_rotation_axis = std::get<4>(initial_data);
        const auto& initial_is_particle_valid = std::get<5>(initial_data);

        const auto timestep_delta = (this->ForceFixedTimeStep || std::get<0>(initial_data) == static_cast<float_t>(0.0))
            ? static_cast<float_t>(this->StreamTimeStep) : std::get<0>(initial_data);

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

        flow_field_module.set_input(*initial_velocities, *initial_global_velocities, initial_get_translation,
            initial_get_rotation_axis, initial_is_particle_valid, seed);
        flow_field_module.set_output(lines);

        flow_field_module.set_parameters(static_cast<tpf::modules::flow_field_aux::method_t>(this->Method),
            static_cast<std::size_t>(this->NumAdvections), timestep_delta);

        flow_field_module.set_callbacks(&call_back_loader);

        flow_field_module.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        auto mesh = tpf::vtk::create_mesh(lines.get_geometry());

        output->SetPoints(mesh.points);
        output->SetLines(mesh.lines);

        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::CELL_DATA, *lines.template get_cell_data_as<std::size_t, 1>("ID (Advection)"));
        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::CELL_DATA, *lines.template get_cell_data_as<std::size_t, 1>("ID (Distribution)"));
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Flow Field' failed."));

        return 0;
    }

    return 1;
}
