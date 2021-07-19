#include "tpf_interface_gradient.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_gradient/tpf_module_interface_gradient.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <exception>

vtkStandardNewMacro(tpf_interface_gradient);

tpf_interface_gradient::tpf_interface_gradient() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_gradient::~tpf_interface_gradient() {}

int tpf_interface_gradient::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_interface_gradient::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::interface_gradient<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_interface_gradient::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_vof = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_vof, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_vof));

        // Create output data
        auto gradients = vof.copy_structure<float_t, 3>("Interface Gradient");

        // Run interface gradient module
        tpf::modules::interface_gradient<float_t> interface_gradient;

        interface_gradient.set_input(vof);
        interface_gradient.set_output(gradients);

        interface_gradient.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->ShallowCopy(in_vof);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, gradients.get_name(), gradients.get_data(), gradients.get_num_components());
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Gradient' failed."));

        return 0;
    }
    
    return 1;
}
