#include "tpf_interface_curvature.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_curvature/tpf_module_interface_curvature.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <stdexcept>
#include <tuple>

vtkStandardNewMacro(tpf_interface_curvature);

tpf_interface_curvature::tpf_interface_curvature() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(3);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_curvature::~tpf_interface_curvature() {}

int tpf_interface_curvature::FillInputPortInformation(int port, vtkInformation* info)
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
    else if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }
    else if (port == 2)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }

    return 0;
}

int tpf_interface_curvature::RequestUpdateExtent(vtkInformation *vtkNotUsed(request), vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::interface_curvature<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_interface_curvature::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));
        const auto gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, in_grid));
        const auto positions = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(2, in_grid));

        // Create output data
        auto curvature = vof.copy_structure<float_t, 1>("Interface Curvature");

        // Run interface gradient module
        tpf::modules::interface_curvature<float_t> interface_curvature;

        interface_curvature.set_input(vof, gradients, positions);
        interface_curvature.set_output(curvature);

        interface_curvature.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->CopyStructure(in_grid);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, curvature.get_name(), curvature.get_data(), curvature.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Curvature' failed."));

        return 0;
    }

    return 1;
}
