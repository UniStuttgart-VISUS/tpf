#include "tpf_temperature.h"

#include "tpf/log/tpf_log.h"

#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <exception>
#include <stdexcept>
#include <string>

vtkStandardNewMacro(tpf_temperature);

tpf_temperature::tpf_temperature()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_temperature::~tpf_temperature() {}

int tpf_temperature::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
        return 1;
    }

    return 0;
}

int tpf_temperature::RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    return 1;
}

int tpf_temperature::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data and information
        auto in_octree = vtkDataSet::GetData(input_vector[0]);

        const auto get_array_name = [this, &in_octree](const int index) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, in_octree);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")."));
        };

        // Get point data
        auto* in_density = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0).c_str()));
        auto* in_density_1 = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1).c_str()));
        auto* in_density_2 = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(2).c_str()));
        auto* in_density_3 = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(3).c_str()));
        auto* in_density_4 = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(4).c_str()));
        auto* in_density_5 = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(5).c_str()));
        auto* in_total_energy = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(6).c_str()));
        auto* in_momentum = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(7).c_str()));
        auto* in_tracer = vtkDoubleArray::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(8).c_str()));
        auto* in_molecular_mass = vtkDoubleArray::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(9).c_str()));
        auto* in_atomic_number = vtkDoubleArray::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(10).c_str()));

        // Sanity checks
        if (in_density->GetNumberOfComponents() != 1 || in_density_1->GetNumberOfComponents() != 1 || in_density_2->GetNumberOfComponents() != 1 ||
            in_density_3->GetNumberOfComponents() != 1 || in_density_4->GetNumberOfComponents() != 1 || in_density_5->GetNumberOfComponents() != 1 ||
            in_total_energy->GetNumberOfComponents() != 1 || in_tracer->GetNumberOfComponents() != 1)
        {
            throw std::runtime_error(__tpf_error_message("Densities and energies must be scalar values."));
        }

        if (in_momentum->GetNumberOfComponents() != 3)
        {
            throw std::runtime_error(__tpf_error_message("Momenta must be three-dimensional vectors."));
        }

        if (in_molecular_mass->GetNumberOfValues() != 5 || in_atomic_number->GetNumberOfValues() != 5)
        {
            throw std::runtime_error(__tpf_error_message("The number of molecular mass and atomic number values must be 5."));
        }

        // Compute temperature: Sagiv Shiber, Simulating White Dwarfs in Octo-Tiger, 2022.
        auto* out_temperature = vtkDoubleArray::New();
        out_temperature->SetName("T");
        out_temperature->SetNumberOfComponents(1);
        out_temperature->SetNumberOfTuples(in_density->GetNumberOfValues());

        const double A = 6.00228 * std::pow(10.0, 22.0);
        const double B = 1.96202 * std::pow(10.0, 6.0);

        const auto molecular_mass_1 = in_molecular_mass->GetValue(0);
        const auto molecular_mass_2 = in_molecular_mass->GetValue(1);
        const auto molecular_mass_3 = in_molecular_mass->GetValue(2);
        const auto molecular_mass_4 = in_molecular_mass->GetValue(3);
        const auto molecular_mass_5 = in_molecular_mass->GetValue(4);

        const auto atomic_number_1 = in_atomic_number->GetValue(0);
        const auto atomic_number_2 = in_atomic_number->GetValue(1);
        const auto atomic_number_3 = in_atomic_number->GetValue(2);
        const auto atomic_number_4 = in_atomic_number->GetValue(3);
        const auto atomic_number_5 = in_atomic_number->GetValue(4);

        const auto mu_1 = molecular_mass_1 / (atomic_number_1 + 1.0);
        const auto mu_2 = molecular_mass_2 / (atomic_number_2 + 1.0);
        const auto mu_3 = molecular_mass_3 / (atomic_number_3 + 1.0);
        const auto mu_4 = molecular_mass_4 / (atomic_number_4 + 1.0);
        const auto mu_5 = molecular_mass_5 / (atomic_number_5 + 1.0);

        const auto hydrogen_mass = 2.01568;
        const auto gamma = 5.0 / 3.0;
        const auto K = 1.380649 * std::pow(10.0, -16.0); // Boltzmann constant

        auto max_E = 0.0;
        auto max_E_deg = 0.0;
        auto max_T = 0.0;
        auto min_T = std::numeric_limits<double>::max();

        for (std::size_t i = 0; i < out_temperature->GetNumberOfValues(); ++i)
        {
            const auto density = in_density->GetValue(i);
            const auto density_1 = in_density_1->GetValue(i);
            const auto density_2 = in_density_2->GetValue(i);
            const auto density_3 = in_density_3->GetValue(i);
            const auto density_4 = in_density_4->GetValue(i);
            const auto density_5 = in_density_5->GetValue(i);
            const auto total_energy = in_total_energy->GetValue(i);
            const auto momentum = in_momentum->GetTuple(i);
            const auto tracer = in_tracer->GetValue(i);

            const auto x = std::pow(density / B, 1.0 / 3.0);
            const auto x_sqr = x * x;
            const auto x_five = std::pow(x, 5.0);
            const auto x_sqr_plus_one_sqrt = std::sqrt(x_sqr + 1);

            const auto P_deg_min = 1.6 * A * x_five;
            const auto P_deg_max = A * (x * (2.0 * x_sqr - 3) * x_sqr_plus_one_sqrt + 3.0 * std::asinh(x));
            const auto P_deg = (x < 0.01) ? P_deg_min : P_deg_max;

            const auto H_deg = (8.0 * A / B) * (x_sqr_plus_one_sqrt - 1.0);

            const auto E_deg_min = 2.4 * A * x_five;
            const auto E_deg_max = H_deg * density - P_deg;
            const auto E_deg = (x < 0.01) ? E_deg_min : E_deg_max;

            const auto E_kin = (0.5 / density) * (momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]);

            const auto E_sub = total_energy - E_kin - E_deg;

            const auto E_th_min = density * std::pow(tracer, gamma);
            const auto E_th_max = E_sub;
            const auto E_th = (E_sub < 0.001 * total_energy) ? E_th_min : E_th_max;

            const auto n = (1.0 / hydrogen_mass) * (
                (density_1 / mu_1) +
                (density_2 / mu_2) +
                (density_3 / mu_3) +
                (density_4 / mu_4) +
                (density_5 / mu_5));

            const auto T = ((gamma - 1.0) * E_th) / (n * K);

            out_temperature->SetValue(i, T);

            // Debug output
            max_E = std::max(max_E, total_energy);
            max_E_deg = std::max(max_E_deg, E_deg);
            min_T = std::min(min_T, T);
            max_T = std::max(max_T, T);
        }

        tpf::log::info_message(__tpf_info_message("max. E_tot: ", max_E, "\nmax. E_deg: ", max_E_deg, "\nT: [", min_T, ", ", max_T, "]\n"));

        // Set output
        auto output = vtkDataSet::GetData(output_vector);
        output->ShallowCopy(in_octree);
        output->GetPointData()->AddArray(out_temperature);
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Temperature (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
