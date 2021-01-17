#include "tpf_flow_field_octree.h"

#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/flow_field/tpf_module_flow_field.h"

#include "tpf/mpi/tpf_mpi.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkDoubleArray.h"
#include "vtkLongLongArray.h"
#define FLOAT_TYPE_ARRAY vtkDoubleArray
#if _WIN32
#define ID_TYPE_ARRAY vtkLongLongArray
#else
#define ID_TYPE_ARRAY vtkLongArray
#endif

#include "vtkObjectFactory.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
    /// Call back class for next time step
    /// <template name="float_t">Floating point type</template>
    /// <template name="point_t">Point type</template>
    template <typename float_t>
    class data_handler : public tpf::modules::flow_field_aux::request_frame_call_back<float_t, tpf::geometry::point<float_t>>
    {
    public:
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="current_timestep">Current time step this plugin was executed in</param>
        /// <param name="timesteps">Available time steps</param>
        /// <param name="request">Original request that will be modified for requesting different time frames</param>
        /// <param name="velocity_alg">Algorithm producing velocity data</param>
        /// <param name="get_array_name">Function to retrieve the array names</param>
        /// <param name="fixed_timestep">Optional fixed timestep</param>
        /// <param name="fixed_frequency">Optional fixed grid rotation</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps, vtkInformation* request,
            vtkAlgorithm* velocity_alg, std::function<std::string(int)> get_array_name,
            std::optional<float_t> fixed_timestep = std::nullopt, std::optional<float_t> fixed_frequency = std::nullopt)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->velocity_alg = velocity_alg;

            this->get_array_name = get_array_name;

            this->fixed_timestep = fixed_timestep;
            this->fixed_frequency = fixed_frequency;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, domain rotation, velocities]</returns>
        virtual std::tuple<float_t, Eigen::Matrix<float_t, 3, 3>,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*> operator()() override
        {
            // Get data
            if (this->timesteps.size() > this->time_offset || this->original_time == this->time_offset)
            {
                vtkIdType num_points;
                vtkPoints* points;

                ID_TYPE_ARRAY* in_paths;
                FLOAT_TYPE_ARRAY* in_x_velocities, * in_y_velocities, * in_z_velocities;

                double x_min, x_max, y_min, y_max, z_min, z_max;
                double time;
                double rotational_frequency;
                double domain_rotation;

                if (tpf::mpi::get_instance().get_rank() == 0)
                {
                    auto in_grid = vtkUnstructuredGrid::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));

                    num_points = in_grid->GetPoints()->GetNumberOfPoints();
                    points = in_grid->GetPoints();

                    in_paths = ID_TYPE_ARRAY::SafeDownCast(in_grid->GetPointData()->GetArray(get_array_name(0).c_str()));
                    in_x_velocities = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetPointData()->GetArray(get_array_name(1).c_str()));
                    in_y_velocities = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetPointData()->GetArray(get_array_name(2).c_str()));
                    in_z_velocities = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetPointData()->GetArray(get_array_name(3).c_str()));

                    x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(4).c_str()))->GetValue(0);
                    x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(5).c_str()))->GetValue(0);
                    y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(6).c_str()))->GetValue(0);
                    y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(7).c_str()))->GetValue(0);
                    z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(8).c_str()))->GetValue(0);
                    z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(9).c_str()))->GetValue(0);

                    time = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(10).c_str()))->GetValue(0);

                    if (this->fixed_frequency)
                    {
                        rotational_frequency = *this->fixed_frequency;
                    }
                    else
                    {
                        rotational_frequency = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(11).c_str()))->GetValue(0);
                    }

                    domain_rotation = FLOAT_TYPE_ARRAY::SafeDownCast(in_grid->GetFieldData()->GetArray(get_array_name(12).c_str()))->GetValue(0);
                }

#ifdef __tpf_use_mpi
                if (tpf::mpi::get_instance().check_mpi_status())
                {
                    tpf::mpi::get_instance().broadcast(num_points, 0);

                    if (tpf::mpi::get_instance().get_rank() != 0)
                    {
                        points = vtkPoints::New();
                        points->SetNumberOfPoints(num_points);

                        in_paths = ID_TYPE_ARRAY::New();
                        in_paths->SetNumberOfComponents(1);
                        in_paths->SetNumberOfTuples(num_points);

                        in_x_velocities = FLOAT_TYPE_ARRAY::New();
                        in_x_velocities->SetNumberOfComponents(1);
                        in_x_velocities->SetNumberOfTuples(num_points);

                        in_y_velocities = FLOAT_TYPE_ARRAY::New();
                        in_y_velocities->SetNumberOfComponents(1);
                        in_y_velocities->SetNumberOfTuples(num_points);

                        in_z_velocities = FLOAT_TYPE_ARRAY::New();
                        in_z_velocities->SetNumberOfComponents(1);
                        in_z_velocities->SetNumberOfTuples(num_points);
                    }

                    MPI_Bcast(points->GetVoidPointer(0), num_points * 3, tpf::mpi::mpi_t<double>::value, 0, tpf::mpi::get_instance().get_comm());

                    MPI_Bcast(in_paths->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

                    MPI_Bcast(in_x_velocities->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
                    MPI_Bcast(in_y_velocities->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
                    MPI_Bcast(in_z_velocities->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

                    tpf::mpi::get_instance().broadcast(x_min, 0);
                    tpf::mpi::get_instance().broadcast(x_max, 0);
                    tpf::mpi::get_instance().broadcast(y_min, 0);
                    tpf::mpi::get_instance().broadcast(y_max, 0);
                    tpf::mpi::get_instance().broadcast(z_min, 0);
                    tpf::mpi::get_instance().broadcast(z_max, 0);

                    tpf::mpi::get_instance().broadcast(time, 0);

                    tpf::mpi::get_instance().broadcast(rotational_frequency, 0);
                    tpf::mpi::get_instance().broadcast(domain_rotation, 0);
                }
#endif

                // Output dataset information
                const auto* vx_range = in_x_velocities->GetValueRange();
                const auto* vy_range = in_y_velocities->GetValueRange();
                const auto* vz_range = in_z_velocities->GetValueRange();

                tpf::log::info_message(__tpf_info_message("Number of points: ", num_points));
                tpf::log::info_message(__tpf_info_message("Domain: [", x_min, ", ", y_min, ", ", z_min, "] x [", x_max, ", ", y_max, ", ", z_max, "]"));
                tpf::log::info_message(__tpf_info_message("Velocities: [", vx_range[0], ", ", vy_range[0], ", ", vz_range[0], "] x [", vx_range[1], ", ", vy_range[1], ", ", vz_range[1], "]"));

                // Create octree with root using the bounding box from the data
                const tpf::geometry::point<float_t> min_point(x_min, y_min, z_min);
                const tpf::geometry::point<float_t> max_point(x_max, y_max, z_max);

                auto bounding_box = std::make_shared<tpf::geometry::cuboid<float_t>>(min_point, max_point);
                auto root = std::make_shared<typename tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>>::node>(std::make_pair(bounding_box, Eigen::Matrix<float_t, 3, 1>()));

                this->octree.set_root(root);

                // Add nodes
                std::vector<std::pair<std::vector<tpf::data::position_t>, Eigen::Matrix<float_t, 3, 1>>> node_information;
                node_information.reserve(num_points);

                for (vtkIdType p = 0; p < num_points; ++p)
                {
                    // Decode path
                    auto path_code = in_paths->GetValue(p);

                    std::vector<tpf::data::position_t> path(static_cast<std::size_t>(std::floor(std::log2(path_code) / std::log2(8))), tpf::data::position_t::INVALID);
                    std::size_t index = 0;

                    while (path_code != 1)
                    {
                        const auto turn = path_code % 8;

                        path[index++] = static_cast<tpf::data::position_t>(((turn & 1) ? 4 : 0) + ((turn & 2) ? 2 : 0) + ((turn & 4) ? 1 : 0));

                        path_code /= 8;
                    }

                    // Add node information
                    node_information.push_back(std::make_pair(path, Eigen::Matrix<float_t, 3, 1>(
                        static_cast<float_t>(in_x_velocities->GetValue(p) - points->GetPoint(p)[1] * rotational_frequency),
                        static_cast<float_t>(in_y_velocities->GetValue(p) + points->GetPoint(p)[0] * rotational_frequency),
                        static_cast<float_t>(in_z_velocities->GetValue(p)))));
                }

                // Delete temporary objects
#ifdef __tpf_use_mpi
                if (tpf::mpi::get_instance().get_rank() != 0)
                {
                    in_paths->Delete();

                    in_x_velocities->Delete();
                    in_y_velocities->Delete();
                    in_z_velocities->Delete();
                }
#endif

                this->octree.insert_nodes(node_information.begin(), node_information.end());

                // Check if next time step is available
                ++this->time_offset;

                if (this->timesteps.size() > this->time_offset)
                {
                    // Request update
                    request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                    this->velocity_alg->Update(request);
                }

                // Calculate time difference
                double next_time;

                if (tpf::mpi::get_instance().get_rank() == 0)
                {
                    auto in_next_grid = vtkUnstructuredGrid::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));
                    next_time = FLOAT_TYPE_ARRAY::SafeDownCast(in_next_grid->GetFieldData()->GetArray(get_array_name(10).c_str()))->GetValue(0);
                }

#ifdef __tpf_use_mpi
                if (tpf::mpi::get_instance().check_mpi_status())
                {
                    tpf::mpi::get_instance().broadcast(next_time, 0);
                }
#endif

                float_t timestep_delta;

                if (this->fixed_timestep)
                {
                    timestep_delta = *this->fixed_timestep;
                }
                else
                {
                    timestep_delta = next_time - time;
                }

                // Calculate rotation matrix for rotation around the z-axis
                Eigen::Matrix<float_t, 3, 3> rotation_matrix;
                rotation_matrix << std::cos(domain_rotation), std::sin(domain_rotation), 0.0,
                                  -std::sin(domain_rotation), std::cos(domain_rotation), 0.0,
                                   0.0,                       0.0,                       0.0;

                return std::make_tuple(timestep_delta, rotation_matrix,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*>(&this->octree));
            }

            throw std::exception();
        }

        /// <summary>
        /// Reset input algorithms to their beginning state
        /// </summary>
        void reset()
        {
            // Request update to original time step
            this->time_offset = this->original_time;

            request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

            this->velocity_alg->Update(request);
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

        /// Velocity algorithm
        vtkAlgorithm* velocity_alg;

        /// Function to retrieve array names
        std::function<std::string(int)> get_array_name;

        /// Fixed timestep
        std::optional<float_t> fixed_timestep;

        /// Fixed rotational frequency
        std::optional<float_t> fixed_frequency;

        /// Data of current timestep
        tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>> octree;
    };
}

vtkStandardNewMacro(tpf_flow_field_octree);

tpf_flow_field_octree::tpf_flow_field_octree()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

tpf_flow_field_octree::~tpf_flow_field_octree() {}

int tpf_flow_field_octree::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }

    return 0;
}

int tpf_flow_field_octree::RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    return 1;
}

int tpf_flow_field_octree::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data and information
        auto in_octree = vtkPointSet::GetData(input_vector[0]);
        auto in_seed = vtkPolyData::GetData(input_vector[1]);

        const auto get_array_name = [this, &in_octree](const int index) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, in_octree);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")"));
        };

        const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_octree).second;
        const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

        data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), get_array_name,
            this->ForceFixedTimeStep ? std::make_optional(static_cast<float_t>(this->StreamTimeStep)) : std::nullopt,
            this->ForceFixedFrequency ? std::make_optional(static_cast<float_t>(this->FrequencyOmega)) : std::nullopt);

        auto initial_data = call_back_loader();
        auto initial_octree = static_cast<tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>>*>(std::get<2>(initial_data));

        const auto timestep_delta = (this->ForceFixedTimeStep || std::get<0>(initial_data) == static_cast<float_t>(0.0))
            ? static_cast<float_t>(this->StreamTimeStep) : std::get<0>(initial_data);

        const auto initial_rotation = std::get<1>(initial_data);

        // Get initial seed
        tpf::data::polydata<float_t> seed = tpf::vtk::get_polydata<float_t>(in_seed);

        // Create and run module for computing integration lines
        using flow_t = tpf::modules::flow_field<float_t, tpf::geometry::point<float_t>>;

        flow_t flow_field;
        tpf::data::polydata<float_t> lines;

        flow_field.set_input(static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>&>(*initial_octree), seed);
        flow_field.set_output(lines);

        flow_field.set_parameters(static_cast<tpf::modules::flow_field_aux::method_t>(this->Method), static_cast<std::size_t>(this->NumAdvections),
            timestep_delta, initial_rotation);

        flow_field.set_callbacks(&call_back_loader);

        flow_field.run();

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
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Flow Field (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
