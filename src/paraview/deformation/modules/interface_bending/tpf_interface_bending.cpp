#include "tpf_interface_bending.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_bending/tpf_module_interface_bending.h"

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

vtkStandardNewMacro(tpf_interface_bending);

tpf_interface_bending::tpf_interface_bending() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_bending::~tpf_interface_bending() {}

int tpf_interface_bending::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_interface_bending::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::interface_bending<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_interface_bending::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto input = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, input));
        const auto gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, input));
        const auto positions = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(2, input));
        const auto velocities = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(3, input));

        // Create output data
        auto min_curvature = vof.copy_structure<float_t, 1>("Bending (minimum)");
        auto max_curvature = vof.copy_structure<float_t, 1>("Bending (maximum)");
        auto absmax_curvature = vof.copy_structure<float_t, 1>("Bending (absolute maximum)");
        auto min_curvature_dir = vof.copy_structure<float_t, 3>("Bending Direction (minimum)");
        auto max_curvature_dir = vof.copy_structure<float_t, 3>("Bending Direction (maximum)");
        auto absmax_curvature_dir = vof.copy_structure<float_t, 3>("Bending Direction (absolute maximum)");

        // Get time step
        const float_t timestep = tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), input);

        // Run interface gradient module
        tpf::modules::interface_bending<float_t> interface_bending;

        interface_bending.set_input(vof, gradients, positions, velocities);
        interface_bending.set_output(min_curvature, max_curvature, absmax_curvature, min_curvature_dir, max_curvature_dir, absmax_curvature_dir);
        interface_bending.set_parameters(timestep);

        interface_bending.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->CopyStructure(input);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, min_curvature.get_name(), min_curvature.get_data(), min_curvature.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, max_curvature.get_name(), max_curvature.get_data(), max_curvature.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, absmax_curvature.get_name(), absmax_curvature.get_data(), absmax_curvature.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, min_curvature_dir.get_name(), min_curvature_dir.get_data(), min_curvature_dir.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, max_curvature_dir.get_name(), max_curvature_dir.get_data(), max_curvature_dir.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, absmax_curvature_dir.get_name(), absmax_curvature_dir.get_data(), absmax_curvature_dir.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Bending' failed."));

        return 0;
    }

    return 1;
}
