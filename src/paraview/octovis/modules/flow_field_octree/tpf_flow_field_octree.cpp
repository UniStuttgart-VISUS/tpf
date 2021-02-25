#include "tpf_flow_field_octree.h"

#include "tpf_modules/flow_field/tpf_module_flow_field.h"

#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/exception/tpf_not_implemented_exception.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_quaternion.h"

#include "tpf/mpi/tpf_mpi.h"

#include "tpf/policies/tpf_interpolatable.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#define FLOAT_TYPE_ARRAY vtkDoubleArray
#ifdef __tpf_no_longlong
#include "vtkLongArray.h"
#define ID_TYPE_ARRAY vtkLongArray
#else
#include "vtkLongLongArray.h"
#define ID_TYPE_ARRAY vtkLongLongArray
#endif

#include "vtkObjectFactory.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

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
        /// <param name="stars_alg">Algorithm producing information about the stars</param>
        /// <param name="get_array_name">Function to retrieve the array names</param>
        /// <param name="locality_method">Method for defining the locality of the stars</param>
        /// <param name="fixed_timestep">Optional fixed timestep</param>
        /// <param name="fixed_frequency">Optional fixed grid angular frequency</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps, vtkInformation* request,
            vtkAlgorithm* velocity_alg, vtkAlgorithm* stars_alg, std::function<std::string(int, int)> get_array_name,
            const tpf_flow_field_octree::locality_method_t locality_method,
            std::optional<float_t> fixed_timestep = std::nullopt, std::optional<float_t> fixed_frequency = std::nullopt)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->velocity_alg = velocity_alg;
            this->stars_alg = stars_alg;

            this->get_array_name = get_array_name;

            this->locality_method = locality_method;

            this->fixed_timestep = fixed_timestep;
            this->fixed_frequency = fixed_frequency;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, velocities, global velocity parts, translation, angular velocity, validity]</returns>
        virtual std::tuple<
            float_t,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<bool(const tpf::geometry::point<float_t>&)>> operator()() override
        {
            // Get data
            if (this->timesteps.size() > this->time_offset || this->original_time == this->time_offset)
            {
                vtkIdType num_points;
                vtkPoints* points;

                ID_TYPE_ARRAY* in_paths, * in_classification;
                FLOAT_TYPE_ARRAY* in_velocities;
                vtkDataArray* in_star_velocities, * in_star_angular_velocity;

                double x_min, x_max, y_min, y_max, z_min, z_max;
                double time;
                double angular_frequency;

                if (tpf::mpi::get_instance().get_rank() == 0)
                {
                    auto in_octree = vtkPointSet::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));
                    auto in_stars = this->stars_alg != nullptr ? vtkPointSet::SafeDownCast(this->stars_alg->GetOutputDataObject(0)) : nullptr;

                    num_points = in_octree->GetPoints()->GetNumberOfPoints();
                    points = in_octree->GetPoints();

                    in_paths = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0, 0).c_str()));

                    in_velocities = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1, 0).c_str()));

                    x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(4, 0).c_str()))->GetValue(0);
                    x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(5, 0).c_str()))->GetValue(0);
                    y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(6, 0).c_str()))->GetValue(0);
                    y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(7, 0).c_str()))->GetValue(0);
                    z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(8, 0).c_str()))->GetValue(0);
                    z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(9, 0).c_str()))->GetValue(0);

                    time = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(10, 0).c_str()))->GetValue(0);

                    // Set angular velocity according to selected method
                    if (this->locality_method == tpf_flow_field_octree::locality_method_t::none)
                    {
                        angular_frequency = 0.0;
                    }
                    else if (this->fixed_frequency)
                    {
                        angular_frequency = *this->fixed_frequency;
                    }
                    else
                    {
                        angular_frequency = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(11, 0).c_str()))->GetValue(0);
                    }

                    this->star_position[0].setZero();
                    this->star_position[1].setZero();

                    // Load data describing translation and rotation of the stars
                    if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                    {
                        if (in_stars == nullptr)
                        {
                            throw std::runtime_error(__tpf_error_message(
                                "Classification and star information necessary for velocity-based and rigid body locality."));
                        }

                        try
                        {
                            in_classification = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(12, 0).c_str()));

                            if (in_classification == nullptr)
                            {
                                throw std::runtime_error(__tpf_error_message("Classification array is null."));
                            }

                            if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                                this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                            {
                                in_star_velocities = in_stars->GetPointData()->GetArray(get_array_name(13, 1).c_str());

                                if (in_star_velocities == nullptr)
                                {
                                    throw std::runtime_error(__tpf_error_message("Classification or star velocity array is null."));
                                }

                                this->star_velocity[0] << in_star_velocities->GetComponent(0, 0), in_star_velocities->GetComponent(0, 1), in_star_velocities->GetComponent(0, 2);
                                this->star_velocity[1] << in_star_velocities->GetComponent(1, 0), in_star_velocities->GetComponent(1, 1), in_star_velocities->GetComponent(1, 2);
                            }

                            if (this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                                this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                            {
                                in_star_angular_velocity = in_stars->GetPointData()->GetArray(get_array_name(14, 1).c_str());

                                if (in_star_angular_velocity == nullptr)
                                {
                                    throw std::runtime_error(__tpf_error_message("Star rotation array is null."));
                                }

                                this->angular_velocity[0] << in_star_angular_velocity->GetComponent(0, 0),
                                    in_star_angular_velocity->GetComponent(0, 1), in_star_angular_velocity->GetComponent(0, 2);
                                this->angular_velocity[1] << in_star_angular_velocity->GetComponent(1, 0),
                                    in_star_angular_velocity->GetComponent(1, 1), in_star_angular_velocity->GetComponent(1, 2);

                                // Get center of mass of the stars
                                std::array<double, 3> tmp_star_pos_1, tmp_star_pos_2;
                                in_stars->GetPoint(0, tmp_star_pos_1.data());
                                in_stars->GetPoint(1, tmp_star_pos_2.data());

                                this->star_position[0] << tmp_star_pos_1[0], tmp_star_pos_1[1], tmp_star_pos_1[2];
                                this->star_position[1] << tmp_star_pos_2[0], tmp_star_pos_2[1], tmp_star_pos_2[2];
                            }
                        }
                        catch (const std::exception& e)
                        {
                            throw std::runtime_error(__tpf_nested_error_message(e.what(),
                                "Classification and star information necessary for velocity-based and rigid body locality."));
                        }
                        catch (...)
                        {
                            throw std::runtime_error(__tpf_error_message(
                                "Classification and star information necessary for velocity-based and rigid body locality."));
                        }
                    }
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

                        in_velocities = FLOAT_TYPE_ARRAY::New();
                        in_velocities->SetNumberOfComponents(3);
                        in_velocities->SetNumberOfTuples(num_points);

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            in_classification = ID_TYPE_ARRAY::New();
                            in_classification->SetNumberOfComponents(1);
                            in_classification->SetNumberOfTuples(num_points);

                            if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                                this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                            {
                                if (vtkFloatArray::SafeDownCast(in_star_velocities) != nullptr)
                                {
                                    in_star_velocities = vtkFloatArray::New();
                                    in_star_velocities->SetNumberOfComponents(3);
                                    in_star_velocities->SetNumberOfTuples(2);
                                }
                                else if (vtkDoubleArray::SafeDownCast(in_star_velocities) != nullptr)
                                {
                                    in_star_velocities = vtkDoubleArray::New();
                                    in_star_velocities->SetNumberOfComponents(3);
                                    in_star_velocities->SetNumberOfTuples(2);
                                }
                                else
                                {
                                    throw std::runtime_error(__tpf_error_message("Invalid star velocity array type."));
                                }
                            }

                            if (this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                                this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                            {
                                if (vtkFloatArray::SafeDownCast(in_star_angular_velocity) != nullptr)
                                {
                                    in_star_angular_velocity = vtkFloatArray::New();
                                    in_star_angular_velocity->SetNumberOfComponents(3);
                                    in_star_angular_velocity->SetNumberOfTuples(2);
                                }
                                else if (vtkDoubleArray::SafeDownCast(in_star_angular_velocity) != nullptr)
                                {
                                    in_star_angular_velocity = vtkDoubleArray::New();
                                    in_star_angular_velocity->SetNumberOfComponents(3);
                                    in_star_angular_velocity->SetNumberOfTuples(2);
                                }
                                else
                                {
                                    throw std::runtime_error(__tpf_error_message("Invalid star rotation array type."));
                                }
                            }
                        }
                    }

                    MPI_Bcast(points->GetVoidPointer(0), num_points * 3, tpf::mpi::mpi_t<double>::value, 0, tpf::mpi::get_instance().get_comm());

                    MPI_Bcast(in_paths->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

                    MPI_Bcast(in_velocities->GetVoidPointer(0), num_points * 3, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

                    tpf::mpi::get_instance().broadcast(x_min, 0);
                    tpf::mpi::get_instance().broadcast(x_max, 0);
                    tpf::mpi::get_instance().broadcast(y_min, 0);
                    tpf::mpi::get_instance().broadcast(y_max, 0);
                    tpf::mpi::get_instance().broadcast(z_min, 0);
                    tpf::mpi::get_instance().broadcast(z_max, 0);

                    tpf::mpi::get_instance().broadcast(time, 0);

                    tpf::mpi::get_instance().broadcast(angular_frequency, 0);

                    if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                    {
                        MPI_Bcast(in_classification->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            MPI_Bcast(in_star_velocities->GetVoidPointer(0), 6,
                                (vtkFloatArray::SafeDownCast(in_star_velocities) != nullptr) ? tpf::mpi::mpi_t<float>::value : tpf::mpi::mpi_t<double>::value,
                                0, tpf::mpi::get_instance().get_comm());
                        }

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            MPI_Bcast(in_star_angular_velocity->GetVoidPointer(0), 6,
                                (vtkFloatArray::SafeDownCast(in_star_angular_velocity) != nullptr) ? tpf::mpi::mpi_t<float>::value : tpf::mpi::mpi_t<double>::value,
                                0, tpf::mpi::get_instance().get_comm());

                            MPI_Bcast(in_star_positions[0].data(), 3, tpf::mpi::mpi_t<float_t>::value, 0, tpf::mpi::get_instance().get_comm());
                            MPI_Bcast(in_star_positions[1].data(), 3, tpf::mpi::mpi_t<float_t>::value, 0, tpf::mpi::get_instance().get_comm());
                        }
                    }
                }
#endif

                // Output dataset information
                std::array<double, 2> vx_range{}, vy_range{}, vz_range{};
                in_velocities->GetValueRange(vx_range.data(), 0);
                in_velocities->GetValueRange(vy_range.data(), 1);
                in_velocities->GetValueRange(vz_range.data(), 2);

                tpf::log::info_message(__tpf_info_message("Number of points: ", num_points));
                tpf::log::info_message(__tpf_info_message("Domain: [", x_min, ", ", y_min, ", ", z_min, "] x [", x_max, ", ", y_max, ", ", z_max, "]"));
                tpf::log::info_message(__tpf_info_message("Velocities: [", vx_range[0], ", ", vy_range[0], ", ", vz_range[0], "] x [", vx_range[1], ", ", vy_range[1], ", ", vz_range[1], "]"));

                // Create octree with root using the bounding box from the data
                const tpf::geometry::point<float_t> min_point(x_min, y_min, z_min);
                const tpf::geometry::point<float_t> max_point(x_max, y_max, z_max);

                auto bounding_box = std::make_shared<tpf::geometry::cuboid<float_t>>(min_point, max_point);
                auto root = std::make_shared<typename tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>>::node>(std::make_pair(bounding_box, Eigen::Matrix<float_t, 3, 1>()));
                auto root_global = std::make_shared<typename tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>>::node>(std::make_pair(bounding_box, Eigen::Matrix<float_t, 3, 1>()));
                auto root_classification = std::make_shared<typename tpf::data::octree<float_t, int>::node>(std::make_pair(bounding_box, 0));

                this->octree.set_root(root);
                this->octree_global.set_root(root_global);

                if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                {
                    this->octree_classification.set_root(root_classification);
                }

                // Add nodes
                std::vector<std::pair<std::vector<tpf::data::position_t>, Eigen::Matrix<float_t, 3, 1>>> node_information, node_information_global;
                std::vector<std::pair<std::vector<tpf::data::position_t>, int>> node_information_classification;
                node_information.reserve(num_points);
                node_information_global.reserve(num_points);

                if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                {
                    node_information_classification.reserve(num_points);
                }

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
                        static_cast<float_t>(in_velocities->GetComponent(p, 0)),
                        static_cast<float_t>(in_velocities->GetComponent(p, 1)),
                        static_cast<float_t>(in_velocities->GetComponent(p, 2)))));

                    std::array<double, 3> point;
                    points->GetPoint(p, point.data());

                    // Compute global velocity part
                    Eigen::Matrix<float_t, 3, 1> tangential_velocity;
                    tangential_velocity.setZero();

                    if ((this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body) &&
                        !this->angular_velocity[in_classification->GetValue(p) - 1].isZero())
                    {
                        const Eigen::Matrix<float_t, 3, 1> position(point[0], point[1], point[2]);
                        const Eigen::Matrix<float_t, 3, 1> relative_position = position - this->star_position[in_classification->GetValue(p) - 1];
                        const Eigen::Matrix<float_t, 3, 1> angular_velocity = this->angular_velocity[in_classification->GetValue(p) - 1];
                        const Eigen::Matrix<float_t, 3, 1> angular_position = relative_position
                            - (relative_position.dot(angular_velocity) / angular_velocity.squaredNorm()) * angular_velocity;

                        tangential_velocity = angular_velocity.cross(angular_position);
                    }

                    switch (this->locality_method)
                    {
                    case tpf_flow_field_octree::locality_method_t::velocity:
                        node_information_global.push_back(std::make_pair(path, this->star_velocity[in_classification->GetValue(p) - 1]));

                        break;
                    case tpf_flow_field_octree::locality_method_t::rotation:
                        node_information_global.push_back(std::make_pair(path, tangential_velocity));

                        break;
                    case tpf_flow_field_octree::locality_method_t::rigid_body:
                        node_information_global.push_back(std::make_pair(path, this->star_velocity[in_classification->GetValue(p) - 1] + tangential_velocity));

                        break;
                    case tpf_flow_field_octree::locality_method_t::none:
                    case tpf_flow_field_octree::locality_method_t::orbit:
                        node_information_global.push_back(std::make_pair(path, Eigen::Matrix<float_t, 3, 1>(
                            static_cast<float_t>(-point[1] * angular_frequency),
                            static_cast<float_t>(point[0] * angular_frequency),
                            static_cast<float_t>(0.0))));
                    }

                    // Add classification information
                    if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                    {
                        node_information_classification.push_back(std::make_pair(path, in_classification->GetValue(p)));
                    }
                }

                // Delete temporary objects
#ifdef __tpf_use_mpi
                if (tpf::mpi::get_instance().get_rank() != 0)
                {
                    in_paths->Delete();

                    in_velocities->Delete();

                    if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                        this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                    {
                        in_classification->Delete();

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            in_star_velocities->Delete();
                        }

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            in_star_angular_velocity->Delete();
                        }
                    }
                }
#endif

                this->octree.insert_nodes(node_information.begin(), node_information.end());
                this->octree_global.insert_nodes(node_information_global.begin(), node_information_global.end());

                if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                {
                    this->octree_classification.insert_nodes(node_information_classification.begin(), node_information_classification.end());
                }

                // Check if next time step is available
                ++this->time_offset;

                if (this->timesteps.size() > this->time_offset)
                {
                    // Request update
                    this->request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                    this->velocity_alg->Update(this->request);
                    this->stars_alg->Update(this->request);
                }

                // Calculate time difference
                double next_time;

                if (tpf::mpi::get_instance().get_rank() == 0)
                {
                    auto in_next_grid = vtkPointSet::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));
                    next_time = FLOAT_TYPE_ARRAY::SafeDownCast(in_next_grid->GetFieldData()->GetArray(get_array_name(10, 0).c_str()))->GetValue(0);
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

                if (timestep_delta == 0.0)
                {
                    timestep_delta = static_cast<float_t>(1.0);
                }

                tpf::log::info_message(__tpf_info_message("Initial time step size: ", timestep_delta));

                // Create functions for locality information
                auto get_star = [this](const tpf::geometry::point<float_t>& position) -> int
                {
                    return this->octree_classification.find_node(position).first->get_value().second;
                };

                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_translation;
                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_angular_velocity;
                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_barycenter;
                std::function<bool(const tpf::geometry::point<float_t>&)> is_valid;

                get_barycenter = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_position[get_star(position) - 1]; };

                switch (this->locality_method)
                {
                case tpf_flow_field_octree::locality_method_t::velocity:
                    get_translation = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_velocity[get_star(position) - 1]; };
                    get_angular_velocity = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::rotation:
                    get_translation = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    get_angular_velocity = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[get_star(position) - 1]; };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::rigid_body:
                    get_translation = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_velocity[get_star(position) - 1]; };
                    get_angular_velocity = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[get_star(position) - 1]; };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::orbit:
                case tpf_flow_field_octree::locality_method_t::none:
                default:
                    this->angular_velocity[0] = this->angular_velocity[1] = Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0) * angular_frequency;

                    get_translation = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    get_angular_velocity = [this](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[0]; };
                    is_valid = [](const tpf::geometry::point<float_t>& position) { return true; };

                    break;
                }

                // Return [Time step delta, velocities, global velocity parts, translation, angular velocity, validity]
                return std::make_tuple(timestep_delta,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*>(&this->octree),
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*>(&this->octree_global),
                    get_translation, get_angular_velocity, get_barycenter, is_valid);
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

            this->velocity_alg->Update(request);
            this->stars_alg->Update(request);
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

        /// Star identification algorithm
        vtkAlgorithm* stars_alg;

        /// Function to retrieve array names
        std::function<std::string(int, int)> get_array_name;

        /// Locality definition for the stars
        typename tpf_flow_field_octree::locality_method_t locality_method;

        /// Fixed timestep
        std::optional<float_t> fixed_timestep;

        /// Fixed rotational frequency
        std::optional<float_t> fixed_frequency;

        /// Data of current timestep
        tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>> octree, octree_global;
        tpf::data::octree<float_t, int> octree_classification;

        /// Position, angular velocity and velocity
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> star_position;
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> angular_velocity;
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> star_velocity;
    };
}

vtkStandardNewMacro(tpf_flow_field_octree);

tpf_flow_field_octree::tpf_flow_field_octree()
{
    this->SetNumberOfInputPorts(3);
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
    if (port == 2)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
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
        auto in_stars = vtkPointSet::GetData(input_vector[2]);

        const auto get_array_name = [this, &in_octree, &in_stars](const int index, const int data_index = 0) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, data_index == 0 ? in_octree : in_stars);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")"));
        };

        const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_octree).second;
        const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

        data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(2, 0),
            get_array_name, static_cast<locality_method_t>(this->LocalityMethod),
            this->ForceFixedTimeStep ? std::make_optional(static_cast<float_t>(this->StreamTimeStep)) : std::nullopt,
            this->ForceFixedFrequency ? std::make_optional(static_cast<float_t>(this->FrequencyOmega)) : std::nullopt);

        // Get initial seed
        tpf::data::polydata<float_t> seed = tpf::vtk::get_polydata<float_t>(in_seed);

        // Create and run module for computing integration lines
        using flow_t = tpf::modules::flow_field<float_t, tpf::geometry::point<float_t>>;

        flow_t flow_field;
        tpf::data::polydata<float_t> lines;

        flow_field.set_input(seed);
        flow_field.set_output(lines);

        flow_field.set_parameters(static_cast<tpf::modules::flow_field_aux::method_t>(this->Method),
            static_cast<std::size_t>(this->NumAdvections));

        flow_field.set_callbacks(&call_back_loader);

        flow_field.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        auto mesh = tpf::vtk::create_mesh(lines.get_geometry());

        output->SetPoints(mesh.points);
        output->SetLines(mesh.lines);

        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::CELL_DATA, *lines.template get_cell_data_as<std::size_t, 1>("ID (Advection)"));
        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::CELL_DATA, *lines.template get_cell_data_as<std::size_t, 1>("ID (Distribution)"));

        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::POINT_DATA, *lines.template get_point_data_as<std::size_t, 1>("ID (Advection)"));
        tpf::vtk::set_data<std::size_t>(output, tpf::data::topology_t::POINT_DATA, *lines.template get_point_data_as<std::size_t, 1>("ID (Distribution)"));
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Flow Field (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
