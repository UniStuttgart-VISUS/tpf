#include "tpf_binary_star.h"

#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_polydata.h"

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
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename float_t>
std::tuple<tpf::data::octree<float_t, float_t>, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*>
load_data(vtkPointSet* in_octree, const std::function<std::string(int)> get_array_name)
{
    // Get input data
    vtkIdType num_points;

    ID_TYPE_ARRAY* in_paths = nullptr;
    FLOAT_TYPE_ARRAY* in_density = nullptr;
    FLOAT_TYPE_ARRAY* in_density_1 = nullptr;
    FLOAT_TYPE_ARRAY* in_density_2 = nullptr;
    FLOAT_TYPE_ARRAY* in_velocity = nullptr;
    FLOAT_TYPE_ARRAY* in_gravitation = nullptr;

    double x_min, x_max, y_min, y_max, z_min, z_max;

    if (tpf::mpi::get_instance().get_rank() == 0)
    {
        num_points = in_octree->GetPoints()->GetNumberOfPoints();

        in_paths = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0).c_str()));

        in_density = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1).c_str()));
        in_density_1 = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(2).c_str()));
        in_density_2 = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(3).c_str()));
        in_velocity = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(4).c_str()));
        in_gravitation = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(5).c_str()));

        x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(6).c_str()))->GetValue(0);
        x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(7).c_str()))->GetValue(0);
        y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(8).c_str()))->GetValue(0);
        y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(9).c_str()))->GetValue(0);
        z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(10).c_str()))->GetValue(0);
        z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(11).c_str()))->GetValue(0);
    }

#ifdef __tpf_use_mpi
    if (tpf::mpi::get_instance().check_mpi_status())
    {
        tpf::mpi::get_instance().broadcast(num_points, 0);

        if (tpf::mpi::get_instance().get_rank() != 0)
        {
            in_paths = ID_TYPE_ARRAY::New();
            in_paths->SetNumberOfComponents(1);
            in_paths->SetNumberOfTuples(num_points);

            in_density = FLOAT_TYPE_ARRAY::New();
            in_density->SetNumberOfComponents(1);
            in_density->SetNumberOfTuples(num_points);

            in_density_1 = FLOAT_TYPE_ARRAY::New();
            in_density_1->SetNumberOfComponents(1);
            in_density_1->SetNumberOfTuples(num_points);

            in_density_2 = FLOAT_TYPE_ARRAY::New();
            in_density_2->SetNumberOfComponents(1);
            in_density_2->SetNumberOfTuples(num_points);

            in_velocity = FLOAT_TYPE_ARRAY::New();
            in_velocity->SetNumberOfComponents(3);
            in_velocity->SetNumberOfTuples(num_points);

            in_gravitation = FLOAT_TYPE_ARRAY::New();
            in_gravitation->SetNumberOfComponents(3);
            in_gravitation->SetNumberOfTuples(num_points);
        }

        MPI_Bcast(in_paths->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

        MPI_Bcast(in_density->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_density_1->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_density_2->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_velocity->GetVoidPointer(0), 3 * num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_gravitation->GetVoidPointer(0), 3 * num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

        tpf::mpi::get_instance().broadcast(x_min, 0);
        tpf::mpi::get_instance().broadcast(x_max, 0);
        tpf::mpi::get_instance().broadcast(y_min, 0);
        tpf::mpi::get_instance().broadcast(y_max, 0);
        tpf::mpi::get_instance().broadcast(z_min, 0);
        tpf::mpi::get_instance().broadcast(z_max, 0);
    }
#endif

    // Create octree with root using the bounding box from the data
    const tpf::geometry::point<float_t> min_point(x_min, y_min, z_min);
    const tpf::geometry::point<float_t> max_point(x_max, y_max, z_max);

    auto bounding_box = std::make_shared<tpf::geometry::cuboid<float_t>>(min_point, max_point);
    auto root = std::make_shared<typename tpf::data::octree<float_t, float_t>::node>(std::make_pair(bounding_box, 0.0f));

    tpf::data::octree<float_t, float_t> octree;
    octree.set_root(root);

    // Add nodes
    std::vector<std::pair<std::vector<tpf::data::position_t>, float_t>> node_information;
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
        node_information.push_back(std::make_pair(path, static_cast<float_t>(in_density->GetValue(p))));
    }

    // Delete temporary objects
#ifdef __tpf_use_mpi
    if (tpf::mpi::get_instance().get_rank() != 0)
    {
        in_paths->Delete();

        in_density->Delete();
        in_density_1->Delete();
        in_density_2->Delete();
        in_velocity->Delete();
        in_gravitation->Delete();
    }
#endif

    octree.insert_nodes(node_information.begin(), node_information.end());

    return std::make_tuple(octree, in_density, in_density_1, in_density_2, in_velocity, in_gravitation);
}

vtkStandardNewMacro(tpf_binary_star);

tpf_binary_star::tpf_binary_star()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_binary_star::~tpf_binary_star() {}

int tpf_binary_star::FillInputPortInformation(int port, vtkInformation* info)
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

    return 0;
}

int tpf_binary_star::RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    return 1;
}

int tpf_binary_star::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data and information
        auto in_octree = vtkPointSet::GetData(input_vector[0]);

        const auto get_array_name = [this, &in_octree](const int index) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, in_octree);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")."));
        };

        // Load data
        const auto data = load_data<float_t>(in_octree, get_array_name);

        const auto octree = std::get<0>(data);

        auto density = std::get<1>(data);
        auto density_1 = std::get<2>(data);
        auto density_2 = std::get<3>(data);
        auto velocities = std::get<4>(data);
        auto gravitation = std::get<5>(data);

        const auto num_points = density->GetNumberOfTuples();

        // Identify donor and accretor cells from the densities
        auto classification = ID_TYPE_ARRAY::New();
        classification->SetName("Classification");
        classification->SetNumberOfComponents(1);
        classification->SetNumberOfTuples(num_points);

        for (vtkIdType p = 0; p < num_points; ++p)
        {
            typename ID_TYPE_ARRAY::ValueType cell_classification;

            // TODO: use density for classification

            classification->SetValue(p, cell_classification);
        }

        // Compute properties for each star or cell
        auto center_of_mass = FLOAT_TYPE_ARRAY::New();
        center_of_mass->SetName("Center of mass");
        center_of_mass->SetNumberOfComponents(3);
        center_of_mass->SetNumberOfTuples(num_points);
        center_of_mass->Fill(0.0);

        auto cluster_velocity = FLOAT_TYPE_ARRAY::New();
        cluster_velocity->SetName("Velocity");
        cluster_velocity->SetNumberOfComponents(3);
        cluster_velocity->SetNumberOfTuples(num_points);
        cluster_velocity->Fill(0.0);

        auto angular_frequency = FLOAT_TYPE_ARRAY::New();
        angular_frequency->SetName("Orbital angular frequency");
        angular_frequency->SetNumberOfComponents(1);
        angular_frequency->SetNumberOfTuples(num_points);
        angular_frequency->Fill(0.0);

        auto acceleration = FLOAT_TYPE_ARRAY::New();
        acceleration->SetName("Acceleration");
        acceleration->SetNumberOfComponents(3);
        acceleration->SetNumberOfTuples(num_points);
        acceleration->Fill(0.0);

        auto classifier = FLOAT_TYPE_ARRAY::New();
        classifier->SetName("Classifier");
        classifier->SetNumberOfComponents(1);
        classifier->SetNumberOfTuples(num_points);
        classifier->Fill(0.0);

        auto roche_lobe = FLOAT_TYPE_ARRAY::New();
        roche_lobe->SetName("Roche lobe radius");
        roche_lobe->SetNumberOfComponents(1);
        roche_lobe->SetNumberOfTuples(num_points);
        roche_lobe->Fill(0.0);

        for (int i = 0; i < 5; ++i)
        {
            // Sum over mass, position and velocity
            std::array<float_t, 2> mass;
            std::array<Eigen::Matrix<float_t, 3, 1>, 2> center;
            std::array<Eigen::Matrix<float_t, 3, 1>, 2> velocity;

            mass[0] = static_cast<float_t>(0.0);
            mass[1] = static_cast<float_t>(0.0);
            center[0].fill(0.0);
            center[1].fill(0.0);
            velocity[0].fill(0.0);
            velocity[1].fill(0.0);

            Eigen::Matrix<double, 3, 1> tmp_position, tmp_velocity;
            Eigen::Matrix<float_t, 3, 1> cell_position, cell_velocity;

            for (vtkIdType p = 0; p < num_points; ++p)
            {
                in_octree->GetPoint(p, tmp_position.data());
                cell_position << static_cast<float_t>(tmp_position[0]),
                    static_cast<float_t>(tmp_position[1]), static_cast<float_t>(tmp_position[2]);

                velocities->GetTuple(p, tmp_velocity.data());
                cell_velocity << static_cast<float_t>(tmp_velocity[0]),
                    static_cast<float_t>(tmp_velocity[1]), static_cast<float_t>(tmp_velocity[2]);

                const auto cell_classification = classification->GetValue(p);
                const auto cell_density = density->GetValue(p);
                const auto cell_volume = octree.find_node(cell_position).first->get_value().first->calculate_volume().get_float_value();

                if (cell_classification != 0)
                {
                    mass[cell_classification - 1] += cell_density * cell_volume;
                    center[cell_classification - 1] += cell_position * cell_density * cell_volume;
                    velocity[cell_classification - 1] += cell_velocity * cell_density * cell_volume;
                }
            }

            // Calculate center of mass and its velocity
            center[0] /= mass[0];
            center[1] /= mass[1];

            velocity[0] /= mass[0];
            velocity[1] /= mass[1];

            // Calculate orbital angular frequency
            const float_t frequency = ((center[0] - center[1]).cross(velocity[0] - velocity[1]) / (center[0] - center[1]).squaredNorm())
                .dot(Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0));

            const auto frequency_squared = frequency * frequency;

            // Calculate Roche lobe radius
            const auto orbital_separation = (center[0] - center[1]).norm();

            const auto ratio_1 = mass[0] / mass[1];
            const auto ratio_2 = mass[1] / mass[0];

            const auto roche_lobe_radius_1 = orbital_separation * static_cast<float_t>(0.49) * std::pow(ratio_1, static_cast<float_t>(2.0 / 3.0)) /
                (static_cast<float_t>(0.6) * std::pow(ratio_1, static_cast<float_t>(2.0 / 3.0))
                    + std::log(static_cast<float_t>(1.0) + std::pow(ratio_1, static_cast<float_t>(1.0 / 3.0))));

            const auto roche_lobe_radius_2 = orbital_separation * static_cast<float_t>(0.49) * std::pow(ratio_2, static_cast<float_t>(2.0 / 3.0)) /
                (static_cast<float_t>(0.6) * std::pow(ratio_2, static_cast<float_t>(2.0 / 3.0))
                    + std::log(static_cast<float_t>(1.0) + std::pow(ratio_2, static_cast<float_t>(1.0 / 3.0))));

            // Iterate over cells, compute the acceleration, and classify it
            Eigen::Matrix<double, 3, 1> tmp_gravitation;
            Eigen::Matrix<float_t, 3, 1> cell_radius, cell_gravitation;

            for (vtkIdType p = 0; p < num_points; ++p)
            {
                in_octree->GetPoint(p, tmp_position.data());
                cell_radius << static_cast<float_t>(tmp_position[0]),
                    static_cast<float_t>(tmp_position[1]), 0.0;
                cell_position << static_cast<float_t>(tmp_position[0]),
                    static_cast<float_t>(tmp_position[1]), static_cast<float_t>(tmp_position[2]);

                velocities->GetTuple(p, tmp_velocity.data());
                cell_velocity << static_cast<float_t>(tmp_velocity[0]),
                    static_cast<float_t>(tmp_velocity[1]), static_cast<float_t>(tmp_velocity[2]);

                gravitation->GetTuple(p, tmp_gravitation.data());
                cell_gravitation << static_cast<float_t>(tmp_gravitation[0]),
                    static_cast<float_t>(tmp_gravitation[1]), static_cast<float_t>(tmp_gravitation[2]);

                const Eigen::Matrix<float_t, 3, 1> cell_acceleration = cell_gravitation + cell_radius * frequency_squared;

                // Classification
                const auto distance_1 = cell_position - center[0];
                const auto distance_2 = cell_position - center[1];

                const float_t cell_classifier_1 = cell_acceleration.dot(distance_1.normalized());
                const float_t cell_classifier_2 = cell_acceleration.dot(distance_2.normalized());

                typename ID_TYPE_ARRAY::ValueType cell_classification;

                if (distance_1.norm() < 0.25 * roche_lobe_radius_1 ||
                    std::min(cell_classifier_1, static_cast<float_t>(0.0)) < std::min(cell_classifier_2, static_cast<float_t>(0.0)))
                {
                    cell_classification = 1;
                }
                else if (distance_2.norm() < 0.25 * roche_lobe_radius_2 ||
                    std::min(cell_classifier_1, static_cast<float_t>(0.0)) > std::min(cell_classifier_2, static_cast<float_t>(0.0)))
                {
                    cell_classification = 2;
                }
                else
                {
                    cell_classification = 0;
                }

                // Store properties
                if (cell_classification != 0)
                {
                    center_of_mass->SetTuple3(p, center[cell_classification - 1][0], center[cell_classification - 1][1], center[cell_classification - 1][2]);
                    cluster_velocity->SetTuple3(p, velocity[cell_classification - 1][0], velocity[cell_classification - 1][1], velocity[cell_classification - 1][2]);
                    angular_frequency->SetValue(p, frequency);
                    acceleration->SetTuple3(p, cell_acceleration[0], cell_acceleration[1], cell_acceleration[2]);
                    classifier->SetValue(p, cell_classification == 1 ? cell_classifier_1 : cell_classifier_2);
                    roche_lobe->SetValue(p, cell_classification == 1 ? roche_lobe_radius_1 : roche_lobe_radius_2);
                }

                classification->SetValue(p, cell_classification);
            }
        }

        // Set output
        auto output = vtkPolyData::GetData(output_vector);
        output->ShallowCopy(in_octree);

        output->GetPointData()->AddArray(center_of_mass);
        output->GetPointData()->AddArray(cluster_velocity);
        output->GetPointData()->AddArray(angular_frequency);
        output->GetPointData()->AddArray(acceleration);
        output->GetPointData()->AddArray(classifier);
        output->GetPointData()->AddArray(roche_lobe);
        output->GetPointData()->AddArray(classification);
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Binary Star (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
