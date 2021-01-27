#include "tpf_binary_star.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"

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
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename float_t>
std::tuple<tpf::data::octree<float_t, float_t>, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*,
    FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*, FLOAT_TYPE_ARRAY*>
load_data(vtkPointSet* in_octree, const std::function<std::string(int)> get_array_name)
{
    // Get input data
    vtkIdType num_points;

    ID_TYPE_ARRAY* in_paths = nullptr;
    FLOAT_TYPE_ARRAY* in_density = nullptr;
    FLOAT_TYPE_ARRAY* in_density_accretor = nullptr;
    FLOAT_TYPE_ARRAY* in_density_donor = nullptr;
    FLOAT_TYPE_ARRAY* in_density_other = nullptr;
    FLOAT_TYPE_ARRAY* in_velocity = nullptr;
    FLOAT_TYPE_ARRAY* in_gravitation = nullptr;

    double x_min, x_max, y_min, y_max, z_min, z_max;

    if (tpf::mpi::get_instance().get_rank() == 0)
    {
        num_points = in_octree->GetPoints()->GetNumberOfPoints();

        in_paths = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0).c_str()));

        in_density = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1).c_str()));
        in_density_accretor = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(2).c_str()));
        in_density_donor = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(3).c_str()));
        in_density_other = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(4).c_str()));
        in_velocity = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(5).c_str()));
        in_gravitation = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(6).c_str()));

        x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(7).c_str()))->GetValue(0);
        x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(8).c_str()))->GetValue(0);
        y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(9).c_str()))->GetValue(0);
        y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(10).c_str()))->GetValue(0);
        z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(11).c_str()))->GetValue(0);
        z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(12).c_str()))->GetValue(0);
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

            in_density_accretor = FLOAT_TYPE_ARRAY::New();
            in_density_accretor->SetNumberOfComponents(1);
            in_density_accretor->SetNumberOfTuples(num_points);

            in_density_donor = FLOAT_TYPE_ARRAY::New();
            in_density_donor->SetNumberOfComponents(1);
            in_density_donor->SetNumberOfTuples(num_points);

            in_density_other = FLOAT_TYPE_ARRAY::New();
            in_density_other->SetNumberOfComponents(1);
            in_density_other->SetNumberOfTuples(num_points);

            in_velocity = FLOAT_TYPE_ARRAY::New();
            in_velocity->SetNumberOfComponents(3);
            in_velocity->SetNumberOfTuples(num_points);

            in_gravitation = FLOAT_TYPE_ARRAY::New();
            in_gravitation->SetNumberOfComponents(3);
            in_gravitation->SetNumberOfTuples(num_points);
        }

        MPI_Bcast(in_paths->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

        MPI_Bcast(in_density->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_density_accretor->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_density_donor->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        MPI_Bcast(in_density_other->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
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
        in_density_accretor->Delete();
        in_density_donor->Delete();
        in_density_other->Delete();
        in_velocity->Delete();
        in_gravitation->Delete();
    }
#endif

    octree.insert_nodes(node_information.begin(), node_information.end());

    return std::make_tuple(octree, in_density, in_density_accretor, in_density_donor, in_density_other, in_velocity, in_gravitation);
}

vtkStandardNewMacro(tpf_binary_star);

tpf_binary_star::tpf_binary_star()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
}

tpf_binary_star::~tpf_binary_star() {}

int tpf_binary_star::ProcessRequest(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Create an output object of the correct type.
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
        return this->RequestDataObject(request, input_vector, output_vector);
    }

    // Generate the data
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
        return this->RequestInformation(request, input_vector, output_vector);
    }

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request, input_vector, output_vector);
    }

    if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
        return this->RequestUpdateExtent(request, input_vector, output_vector);
    }

    return this->Superclass::ProcessRequest(request, input_vector, output_vector);
}

int tpf_binary_star::RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    // Octree
    {
        auto output = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));

        if (!output)
        {
            output = vtkPolyData::New();
            output_vector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), output);
            this->GetOutputPortInformation(0)->Set(vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
        }
    }

    // Stars
    {
        auto output = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));

        if (!output)
        {
            output = vtkPolyData::New();
            output_vector->GetInformationObject(1)->Set(vtkDataObject::DATA_OBJECT(), output);
            this->GetOutputPortInformation(1)->Set(vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
        }
    }

    return 1;
}

int tpf_binary_star::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }

    return 0;
}

int tpf_binary_star::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }
    else if (port == 1)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }

    return 1;
}

int tpf_binary_star::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int tpf_binary_star::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int tpf_binary_star::RequestData(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
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

        // Create stars poly data object
        tpf::data::polydata<float_t> stars;

        stars.insert(std::make_shared<tpf::geometry::point<float_t>>(0.0, 0.0, 0.0));
        stars.insert(std::make_shared<tpf::geometry::point<float_t>>(0.0, 0.0, 0.0));

        stars.add(std::make_shared<tpf::data::array<float_t, 3>>("Velocity", 2), tpf::data::topology_t::POINT_DATA);
        stars.add(std::make_shared<tpf::data::array<float_t, 1>>("Mass", 2), tpf::data::topology_t::POINT_DATA);
        stars.add(std::make_shared<tpf::data::array<float_t, 1>>("Orbital frequency", 2), tpf::data::topology_t::POINT_DATA);
        stars.add(std::make_shared<tpf::data::array<float_t, 1>>("Roche lobe radius", 2), tpf::data::topology_t::POINT_DATA);

        // Load data
        const auto data = load_data<float_t>(in_octree, get_array_name);

        const auto octree = std::get<0>(data);

        auto density = std::get<1>(data);
        auto density_accretor = std::get<2>(data);
        auto density_donor = std::get<3>(data);
        auto density_other = std::get<4>(data);
        auto velocities = std::get<5>(data);
        auto gravitation = std::get<6>(data);

        const auto num_points = density->GetNumberOfTuples();

        // "Enum" for classification
        const typename ID_TYPE_ARRAY::ValueType OTHER = 0;
        const typename ID_TYPE_ARRAY::ValueType ACCRETOR = 1;
        const typename ID_TYPE_ARRAY::ValueType DONOR = 2;

        const typename ID_TYPE_ARRAY::ValueType ACCRETOR_INDEX = 0;
        const typename ID_TYPE_ARRAY::ValueType DONOR_INDEX = 1;

        // Identify donor and accretor cells from the densities
        auto classification = vtkSmartPointer<ID_TYPE_ARRAY>::New();
        classification->SetName("Classification");
        classification->SetNumberOfComponents(1);
        classification->SetNumberOfTuples(num_points);

        for (vtkIdType p = 0; p < num_points; ++p)
        {
            typename ID_TYPE_ARRAY::ValueType cell_classification;

            const auto cell_density = density->GetValue(p);
            const auto cell_density_accretor = density_accretor->GetValue(p);
            const auto cell_density_donor = density_donor->GetValue(p);
            const auto cell_density_other = density_other->GetValue(p);

            if (cell_density < this->DensityCutoff ||
                (cell_density_other > cell_density_accretor && cell_density_other > cell_density_donor))
            {
                cell_classification = OTHER;
            }
            else if (cell_density_accretor > cell_density_donor)
            {
                cell_classification = ACCRETOR;
            }
            else
            {
                cell_classification = DONOR;
            }

            classification->SetValue(p, cell_classification);
        }

        // Compute properties for each star or cell
        auto center_of_mass = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        center_of_mass->SetName("Center of mass");
        center_of_mass->SetNumberOfComponents(3);
        center_of_mass->SetNumberOfTuples(num_points);
        center_of_mass->Fill(0.0);

        auto cluster_velocity = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        cluster_velocity->SetName("Velocity");
        cluster_velocity->SetNumberOfComponents(3);
        cluster_velocity->SetNumberOfTuples(num_points);
        cluster_velocity->Fill(0.0);

        auto angular_frequency = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        angular_frequency->SetName("Orbital angular frequency");
        angular_frequency->SetNumberOfComponents(1);
        angular_frequency->SetNumberOfTuples(num_points);
        angular_frequency->Fill(0.0);

        auto acceleration = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        acceleration->SetName("Acceleration");
        acceleration->SetNumberOfComponents(3);
        acceleration->SetNumberOfTuples(num_points);
        acceleration->Fill(0.0);

        auto classifier = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        classifier->SetName("Classifier");
        classifier->SetNumberOfComponents(1);
        classifier->SetNumberOfTuples(num_points);
        classifier->Fill(0.0);

        auto roche_lobe = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        roche_lobe->SetName("Roche lobe radius");
        roche_lobe->SetNumberOfComponents(1);
        roche_lobe->SetNumberOfTuples(num_points);
        roche_lobe->Fill(0.0);

        const int value_shift = static_cast<int>(std::floor(std::log(octree.get_root()->get_value().first->calculate_volume().get_float_value())));
        const float_t value_multiplier = static_cast<float_t>(std::pow(2.0, -value_shift));

        for (int i = 0; i <= this->NumIterations; ++i)
        {
            // Sum over mass, position and velocity
            std::array<float_t, 2> mass;
            std::array<Eigen::Matrix<float_t, 3, 1>, 2> center;
            std::array<Eigen::Matrix<float_t, 3, 1>, 2> velocity;

            mass[ACCRETOR_INDEX] = static_cast<float_t>(0.0);
            mass[DONOR_INDEX] = static_cast<float_t>(0.0);
            center[ACCRETOR_INDEX].fill(0.0);
            center[DONOR_INDEX].fill(0.0);
            velocity[ACCRETOR_INDEX].fill(0.0);
            velocity[DONOR_INDEX].fill(0.0);

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

                if (cell_classification != OTHER)
                {
                    // Values are too big; make them smaller
                    const auto cell_mass = (cell_density * cell_volume) * value_multiplier;

                    mass[cell_classification - 1] += cell_mass;
                    center[cell_classification - 1] += cell_position * cell_mass;
                    velocity[cell_classification - 1] += cell_velocity * cell_mass;
                }
            }

            // Calculate center of mass and its velocity
            center[ACCRETOR_INDEX] /= mass[ACCRETOR_INDEX];
            center[DONOR_INDEX] /= mass[DONOR_INDEX];

            velocity[ACCRETOR_INDEX] /= mass[ACCRETOR_INDEX];
            velocity[DONOR_INDEX] /= mass[DONOR_INDEX];

            mass[ACCRETOR_INDEX] /= value_multiplier;
            mass[DONOR_INDEX] /= value_multiplier;

            // Calculate orbital angular frequency
            const float_t frequency = ((center[ACCRETOR_INDEX] - center[DONOR_INDEX])
                .cross(velocity[ACCRETOR_INDEX] - velocity[DONOR_INDEX]) / (center[ACCRETOR_INDEX] - center[DONOR_INDEX]).squaredNorm())
                .dot(Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0));

            const auto frequency_squared = frequency * frequency;

            // Calculate Roche lobe radius
            const auto orbital_separation = (center[ACCRETOR_INDEX] - center[DONOR_INDEX]).norm();

            const auto ratio_accretor = mass[ACCRETOR_INDEX] / mass[DONOR_INDEX];
            const auto ratio_donor = mass[DONOR_INDEX] / mass[ACCRETOR_INDEX];

            const auto roche_lobe_radius_accretor = orbital_separation * static_cast<float_t>(0.49) * std::pow(ratio_accretor, static_cast<float_t>(2.0 / 3.0)) /
                (static_cast<float_t>(0.6) * std::pow(ratio_accretor, static_cast<float_t>(2.0 / 3.0))
                    + std::log(static_cast<float_t>(1.0) + std::pow(ratio_accretor, static_cast<float_t>(1.0 / 3.0))));

            const auto roche_lobe_radius_donor = orbital_separation * static_cast<float_t>(0.49) * std::pow(ratio_donor, static_cast<float_t>(2.0 / 3.0)) /
                (static_cast<float_t>(0.6) * std::pow(ratio_donor, static_cast<float_t>(2.0 / 3.0))
                    + std::log(static_cast<float_t>(1.0) + std::pow(ratio_donor, static_cast<float_t>(1.0 / 3.0))));

            if (i == this->NumIterations)
            {
#if __tpf_detailed
                // Output information
                tpf::log::info_message(__tpf_info_message("Accretor"));
                tpf::log::info_message(__tpf_info_message("  mass              ", mass[ACCRETOR_INDEX]));
                tpf::log::info_message(__tpf_info_message("  location          [", center[ACCRETOR_INDEX][0], ", ", center[ACCRETOR_INDEX][1], ", ", center[ACCRETOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  velocity          [", velocity[ACCRETOR_INDEX][0], ", ", velocity[ACCRETOR_INDEX][1], ", ", velocity[ACCRETOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  angular frequency ", frequency));
                tpf::log::info_message(__tpf_info_message("  Roche lobe radius ", roche_lobe_radius_accretor));

                tpf::log::info_message(__tpf_info_message("Donor"));
                tpf::log::info_message(__tpf_info_message("  mass              ", mass[DONOR_INDEX]));
                tpf::log::info_message(__tpf_info_message("  location          [", center[DONOR_INDEX][0], ", ", center[DONOR_INDEX][1], ", ", center[DONOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  velocity          [", velocity[DONOR_INDEX][0], ", ", velocity[DONOR_INDEX][1], ", ", velocity[DONOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  angular frequency ", frequency));
                tpf::log::info_message(__tpf_info_message("  Roche lobe radius ", roche_lobe_radius_donor));
#endif

                // Set stars
                *std::static_pointer_cast<tpf::geometry::point<float_t>>(stars.get_object(ACCRETOR_INDEX)) = center[ACCRETOR_INDEX];
                *std::static_pointer_cast<tpf::geometry::point<float_t>>(stars.get_object(DONOR_INDEX)) = center[DONOR_INDEX];

                stars.template get_point_data_as<float_t, 3, 1>("Velocity")->at(ACCRETOR_INDEX) = velocity[ACCRETOR_INDEX];
                stars.template get_point_data_as<float_t, 3, 1>("Velocity")->at(DONOR_INDEX) = velocity[DONOR_INDEX];

                stars.template get_point_data_as<float_t, 1, 1>("Mass")->at(ACCRETOR_INDEX) = mass[ACCRETOR_INDEX];
                stars.template get_point_data_as<float_t, 1, 1>("Mass")->at(DONOR_INDEX) = mass[DONOR_INDEX];

                stars.template get_point_data_as<float_t, 1, 1>("Orbital frequency")->at(ACCRETOR_INDEX) = frequency;
                stars.template get_point_data_as<float_t, 1, 1>("Orbital frequency")->at(DONOR_INDEX) = frequency;

                stars.template get_point_data_as<float_t, 1, 1>("Roche lobe radius")->at(ACCRETOR_INDEX) = roche_lobe_radius_accretor;
                stars.template get_point_data_as<float_t, 1, 1>("Roche lobe radius")->at(DONOR_INDEX) = roche_lobe_radius_donor;

                break;
            }

            // Iterate over cells, compute the acceleration, and classify it
            Eigen::Matrix<double, 3, 1> tmp_gravitation;
            Eigen::Matrix<float_t, 3, 1> cell_radius, cell_gravitation, cell_acceleration;
            float_t cell_classifier_accretor, cell_classifier_donor;

            for (vtkIdType p = 0; p < num_points; ++p)
            {
                typename ID_TYPE_ARRAY::ValueType cell_classification = OTHER;

                if (density->GetValue(p) > this->DensityCutoff)
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

                    cell_acceleration = cell_gravitation + cell_radius * frequency_squared;

                    // Classification
                    const auto distance_accretor = cell_position - center[ACCRETOR_INDEX];
                    const auto distance_donor = cell_position - center[DONOR_INDEX];

                    cell_classifier_accretor = cell_acceleration.dot(distance_accretor.normalized());
                    cell_classifier_donor = cell_acceleration.dot(distance_donor.normalized());

                    if (distance_accretor.norm() < 0.25 * roche_lobe_radius_accretor ||
                        std::min(cell_classifier_accretor, static_cast<float_t>(0.0)) < std::min(cell_classifier_donor, static_cast<float_t>(0.0)))
                    {
                        cell_classification = ACCRETOR;
                    }
                    else if (distance_donor.norm() < 0.25 * roche_lobe_radius_donor ||
                        std::min(cell_classifier_accretor, static_cast<float_t>(0.0)) > std::min(cell_classifier_donor, static_cast<float_t>(0.0)))
                    {
                        cell_classification = DONOR;
                    }
                }

                // Store properties
                if (cell_classification != OTHER)
                {
                    center_of_mass->SetTuple3(p, center[cell_classification - 1][0], center[cell_classification - 1][1], center[cell_classification - 1][2]);
                    cluster_velocity->SetTuple3(p, velocity[cell_classification - 1][0], velocity[cell_classification - 1][1], velocity[cell_classification - 1][2]);
                    angular_frequency->SetValue(p, frequency);
                    acceleration->SetTuple3(p, cell_acceleration[0], cell_acceleration[1], cell_acceleration[2]);
                    classifier->SetValue(p, cell_classification == 1 ? cell_classifier_accretor : cell_classifier_donor);
                    roche_lobe->SetValue(p, cell_classification == 1 ? roche_lobe_radius_accretor : roche_lobe_radius_donor);
                }

                classification->SetValue(p, cell_classification);
            }
        }

        // Create output array
        auto omega = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        omega->SetName("Orbital frequency");
        omega->SetNumberOfComponents(1);
        omega->SetNumberOfTuples(1);
        omega->SetValue(0, stars.template get_point_data_as<float_t, 1, 1>("Orbital frequency")->at(ACCRETOR_INDEX));

        // Set output octree
        auto output_octree = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
        output_octree->ShallowCopy(in_octree);

        output_octree->GetPointData()->AddArray(center_of_mass);
        output_octree->GetPointData()->AddArray(cluster_velocity);
        output_octree->GetPointData()->AddArray(angular_frequency);
        output_octree->GetPointData()->AddArray(acceleration);
        output_octree->GetPointData()->AddArray(classifier);
        output_octree->GetPointData()->AddArray(roche_lobe);
        output_octree->GetPointData()->AddArray(classification);

        output_octree->GetFieldData()->AddArray(omega);

        // Set output stars
        auto output_stars = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));

        tpf::vtk::set_polydata(output_stars, stars,
            tpf::data::data_information<float_t, 3>{ "Velocity", tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ "Mass", tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ "Orbital frequency", tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ "Roche lobe radius", tpf::data::topology_t::POINT_DATA });

        output_stars->GetFieldData()->AddArray(omega);
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Binary Star (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
