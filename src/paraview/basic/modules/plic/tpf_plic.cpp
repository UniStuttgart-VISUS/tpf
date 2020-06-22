#include "tpf_plic.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_gradient/tpf_module_interface_gradient.h"
#include "tpf_modules/plic/tpf_module_plic.h"

#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <stdexcept>
#include <tuple>

vtkStandardNewMacro(tpf_plic);

tpf_plic::tpf_plic() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_plic::~tpf_plic() {}

int tpf_plic::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_plic::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::plic<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_plic::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, this->GetInputArrayToProcess(0, in_grid));

        // Create output data
        auto plic_interface = tpf::data::polydata<float_t>();

        // Run interface gradient module to calculate the gradient if necessary
        const bool has_gradients = this->GetInputArrayToProcess(1, in_grid) != nullptr;

        auto gradients = vof.template copy_structure<float_t, 3>("Interface Gradient");

        if (!has_gradients)
        {
            tpf::modules::interface_gradient<float_t> interface_gradient_module;

            interface_gradient_module.set_input(vof);
            interface_gradient_module.set_output(gradients);
            interface_gradient_module.run();
        }
        else
        {
            gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, this->GetInputArrayToProcess(1, in_grid));
        }

        // Run PLIC module
        tpf::modules::plic<float_t> plic_module;

        plic_module.set_input(vof, gradients);
        plic_module.set_output(plic_interface);
        plic_module.set_parameters(this->NumIterations, static_cast<float_t>(this->Perturbation));
        plic_module.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

#ifndef __tpf_debug
        tpf::vtk::set_polydata(output, plic_interface);
#else
        tpf::vtk::set_polydata(output, plic_interface, tpf::data::data_information<float_t, 1>{ tpf::module_name::plic::error, tpf::data::topology_t::CELL_DATA });
#endif
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'PLIC' failed."));

        return 0;
    }

    return 1;
}
