#include "tpf_binary_star.h"

#include "tpf_modules/droplets/tpf_cluster.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_position.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkDoubleArray.h"
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

vtkStandardNewMacro(tpf_binary_star);

tpf_binary_star::tpf_binary_star()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
}

tpf_binary_star::~tpf_binary_star() { }

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
        stars.add(std::make_shared<tpf::data::array<float_t, 3>>("Axis of rotation", 2), tpf::data::topology_t::POINT_DATA);
        stars.add(std::make_shared<tpf::data::array<float_t, 1>>("Orbital frequency", 2), tpf::data::topology_t::POINT_DATA);
        stars.add(std::make_shared<tpf::data::array<float_t, 1>>("Roche lobe radius", 2), tpf::data::topology_t::POINT_DATA);

        // Get input data
        vtkIdType num_points;

        std::shared_ptr<FLOAT_TYPE_ARRAY> cell_size = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> density = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> density_accretor = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> density_donor = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> density_other = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> velocities = nullptr;
        std::shared_ptr<FLOAT_TYPE_ARRAY> gravitation = nullptr;

        double x_min, x_max, y_min, y_max, z_min, z_max;

        if (tpf::mpi::get_instance().get_rank() == 0)
        {
            auto deleter = [](FLOAT_TYPE_ARRAY*) {};

            num_points = in_octree->GetPoints()->GetNumberOfPoints();

            cell_size = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0).c_str())), deleter);

            density = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1).c_str())), deleter);
            density_accretor = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(2).c_str())), deleter);
            density_donor = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(3).c_str())), deleter);
            density_other = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(4).c_str())), deleter);
            velocities = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(5).c_str())), deleter);
            gravitation = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(6).c_str())), deleter);

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
                auto deleter = [](FLOAT_TYPE_ARRAY* ptr) { ptr->Delete(); };

                density = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                density->SetNumberOfComponents(1);
                density->SetNumberOfTuples(num_points);

                density_accretor = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                density_accretor->SetNumberOfComponents(1);
                density_accretor->SetNumberOfTuples(num_points);

                density_donor = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                density_donor->SetNumberOfComponents(1);
                density_donor->SetNumberOfTuples(num_points);

                density_other = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                density_other->SetNumberOfComponents(1);
                density_other->SetNumberOfTuples(num_points);

                velocities = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                velocities->SetNumberOfComponents(3);
                velocities->SetNumberOfTuples(num_points);

                gravitation = std::shared_ptr<FLOAT_TYPE_ARRAY>(FLOAT_TYPE_ARRAY::New(), deleter);
                gravitation->SetNumberOfComponents(3);
                gravitation->SetNumberOfTuples(num_points);
            }

            MPI_Bcast(density->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(density_accretor->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(density_donor->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(density_other->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(velocities->GetVoidPointer(0), 3 * num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(gravitation->GetVoidPointer(0), 3 * num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
            MPI_Bcast(cell_size->GetVoidPointer(0), 3 * num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

            tpf::mpi::get_instance().broadcast(x_min, 0);
            tpf::mpi::get_instance().broadcast(x_max, 0);
            tpf::mpi::get_instance().broadcast(y_min, 0);
            tpf::mpi::get_instance().broadcast(y_max, 0);
            tpf::mpi::get_instance().broadcast(z_min, 0);
            tpf::mpi::get_instance().broadcast(z_max, 0);
        }
#endif

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
#ifdef __tpf_debug
        auto center_of_mass = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        center_of_mass->SetName("Center of mass");
        center_of_mass->SetNumberOfComponents(3);
        center_of_mass->SetNumberOfTuples(num_points);
        center_of_mass->Fill(0.0);

        auto cluster_velocity = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        cluster_velocity->SetName("Star velocity");
        cluster_velocity->SetNumberOfComponents(3);
        cluster_velocity->SetNumberOfTuples(num_points);
        cluster_velocity->Fill(0.0);

        auto cluster_rotation = vtkSmartPointer<FLOAT_TYPE_ARRAY>::New();
        cluster_rotation->SetName("Axis of rotation");
        cluster_rotation->SetNumberOfComponents(3);
        cluster_rotation->SetNumberOfTuples(num_points);
        cluster_rotation->Fill(0.0);

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
#endif

        for (int i = 0; i <= this->NumIterations; ++i)
        {
            // Sum over mass, position and velocity
            std::array<tpf::modules::droplets_aux::cluster<float_t>, 2> cluster;

            Eigen::Matrix<double, 3, 1> tmp_position, tmp_velocity;
            Eigen::Matrix<float_t, 3, 1> cell_position, cell_velocity;

            for (vtkIdType p = 0; p < num_points; ++p)
            {
                const auto cell_classification = classification->GetValue(p);

                if (cell_classification != OTHER)
                {
                    in_octree->GetPoint(p, tmp_position.data());
                    cell_position << static_cast<float_t>(tmp_position[0]),
                        static_cast<float_t>(tmp_position[1]), static_cast<float_t>(tmp_position[2]);

                    velocities->GetTuple(p, tmp_velocity.data());
                    cell_velocity << static_cast<float_t>(tmp_velocity[0]),
                        static_cast<float_t>(tmp_velocity[1]), static_cast<float_t>(tmp_velocity[2]);

                    const auto cell_density = density->GetValue(p);
                    const auto cell_volume = std::pow(2.0, -9.0) * cell_size->GetComponent(p, 0) * cell_size->GetComponent(p, 1) * cell_size->GetComponent(p, 2);

                    // Values are too big; make them smaller
                    const auto cell_mass = (cell_density * cell_volume);

                    cluster[cell_classification - 1].add(tpf::data::coords3_t(), cell_mass, cell_position, cell_velocity);
                }
            }

            // Calculate center of mass and its velocity
            const std::array<float_t, 2> mass{ cluster[0].get_volume(), cluster[1].get_volume() };
            const std::array<Eigen::Matrix<float_t, 3, 1>, 2> center{ cluster[0].get_position(), cluster[1].get_position() };
            const std::array<Eigen::Matrix<float_t, 3, 1>, 2> velocity{ cluster[0].get_velocity().first, cluster[1].get_velocity().first };

            // Compute axis of rotation
            const std::array<Eigen::Matrix<float_t, 3, 1>, 2> rotation_axis{
                std::get<0>(cluster[0].get_rotation_mechanics()), std::get<0>(cluster[1].get_rotation_mechanics()) };

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
#ifdef __tpf_detailed
                // Output information
                tpf::log::info_message(__tpf_info_message("Accretor"));
                tpf::log::info_message(__tpf_info_message("  mass              ", mass[ACCRETOR_INDEX]));
                tpf::log::info_message(__tpf_info_message("  location          [", center[ACCRETOR_INDEX][0], ", ", center[ACCRETOR_INDEX][1], ", ", center[ACCRETOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  velocity          [", velocity[ACCRETOR_INDEX][0], ", ", velocity[ACCRETOR_INDEX][1], ", ", velocity[ACCRETOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  axis of rotation  [", rotation_axis[ACCRETOR_INDEX][0], ", ", rotation_axis[ACCRETOR_INDEX][1], ", ", rotation_axis[ACCRETOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  angular frequency ", frequency));
                tpf::log::info_message(__tpf_info_message("  Roche lobe radius ", roche_lobe_radius_accretor));

                tpf::log::info_message(__tpf_info_message("Donor"));
                tpf::log::info_message(__tpf_info_message("  mass              ", mass[DONOR_INDEX]));
                tpf::log::info_message(__tpf_info_message("  location          [", center[DONOR_INDEX][0], ", ", center[DONOR_INDEX][1], ", ", center[DONOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  velocity          [", velocity[DONOR_INDEX][0], ", ", velocity[DONOR_INDEX][1], ", ", velocity[DONOR_INDEX][2], "]"));
                tpf::log::info_message(__tpf_info_message("  axis of rotation  [", rotation_axis[DONOR_INDEX][0], ", ", rotation_axis[DONOR_INDEX][1], ", ", rotation_axis[DONOR_INDEX][2], "]"));
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

                stars.template get_point_data_as<float_t, 3, 1>("Axis of rotation")->at(ACCRETOR_INDEX) = rotation_axis[ACCRETOR_INDEX];
                stars.template get_point_data_as<float_t, 3, 1>("Axis of rotation")->at(DONOR_INDEX) = rotation_axis[DONOR_INDEX];

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
#ifdef __tpf_debug
                if (cell_classification != OTHER)
                {
                    center_of_mass->SetTuple3(p, center[cell_classification - 1][0], center[cell_classification - 1][1], center[cell_classification - 1][2]);
                    cluster_velocity->SetTuple3(p, velocity[cell_classification - 1][0], velocity[cell_classification - 1][1], velocity[cell_classification - 1][2]);
                    cluster_rotation->SetTuple3(p, rotation_axis[cell_classification - 1][0], rotation_axis[cell_classification - 1][1], rotation_axis[cell_classification - 1][2]);
                    angular_frequency->SetValue(p, frequency);
                    acceleration->SetTuple3(p, cell_acceleration[0], cell_acceleration[1], cell_acceleration[2]);
                    classifier->SetValue(p, cell_classification == 1 ? cell_classifier_accretor : cell_classifier_donor);
                    roche_lobe->SetValue(p, cell_classification == 1 ? roche_lobe_radius_accretor : roche_lobe_radius_donor);
                }
#endif

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

#ifdef __tpf_debug
        output_octree->GetPointData()->AddArray(center_of_mass);
        output_octree->GetPointData()->AddArray(cluster_velocity);
        output_octree->GetPointData()->AddArray(cluster_rotation);
        output_octree->GetPointData()->AddArray(angular_frequency);
        output_octree->GetPointData()->AddArray(acceleration);
        output_octree->GetPointData()->AddArray(classifier);
        output_octree->GetPointData()->AddArray(roche_lobe);
#endif

        output_octree->GetPointData()->AddArray(classification);

        output_octree->GetFieldData()->AddArray(omega);

        // Set output stars
        auto output_stars = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));

        tpf::vtk::set_polydata(output_stars, stars,
            tpf::data::data_information<float_t, 3>{ "Velocity", tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ "Mass", tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 3>{ "Axis of rotation", tpf::data::topology_t::POINT_DATA },
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
