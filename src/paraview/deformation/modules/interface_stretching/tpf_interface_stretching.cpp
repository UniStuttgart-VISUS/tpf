#include "tpf_interface_stretching.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_stretching/tpf_module_interface_stretching.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <stdexcept>

vtkStandardNewMacro(tpf_interface_stretching);

tpf_interface_stretching::tpf_interface_stretching() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_stretching::~tpf_interface_stretching() {}

int tpf_interface_stretching::FillInputPortInformation(int port, vtkInformation* info)
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

    return 0;
}

int tpf_interface_stretching::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::interface_stretching<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_interface_stretching::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto input = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, input));
        const auto gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, input));
        const auto positions = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(2, input));
        const auto velocities = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(3, input));

        // Get time step
        const float_t timestep = tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), input);

        // Create output data
        auto stretching = vof.copy_structure<float_t, 1>("Stretching (area)");
        auto stretching_min_dir = vof.copy_structure<float_t, 3>("Stretching Direction (minimum)");
        auto stretching_max_dir = vof.copy_structure<float_t, 3>("Stretching Direction (maximum)");
        auto stretching_largest_dir = vof.copy_structure<float_t, 3>("Stretching Direction (largest)");

        // Run interface gradient module
        tpf::modules::interface_stretching<float_t> interface_stretching;

        interface_stretching.set_input(vof, gradients, positions, velocities);
        interface_stretching.set_output(stretching, stretching_max_dir, stretching_min_dir, stretching_largest_dir);
        interface_stretching.set_parameters(timestep);

        interface_stretching.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->ShallowCopy(input);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretching.get_name(), stretching.get_data(), stretching.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretching_max_dir.get_name(), stretching_max_dir.get_data(), stretching_max_dir.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretching_min_dir.get_name(), stretching_min_dir.get_data(), stretching_min_dir.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretching_largest_dir.get_name(), stretching_largest_dir.get_data(), stretching_largest_dir.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Stretching' failed."));

        return 0;
    }

    return 1;
}
