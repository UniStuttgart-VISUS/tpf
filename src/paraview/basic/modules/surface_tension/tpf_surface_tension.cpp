#include "tpf_surface_tension.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/surface_tension/tpf_module_surface_tension.h"

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

vtkStandardNewMacro(tpf_surface_tension);

tpf_surface_tension::tpf_surface_tension() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(3);
    this->SetNumberOfOutputPorts(1);
}

tpf_surface_tension::~tpf_surface_tension() {}

int tpf_surface_tension::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_surface_tension::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::surface_tension<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_surface_tension::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);

        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));
        const auto gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, in_grid));
        const auto curvature = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(2, in_grid));

        // Create output data
        auto surface_tension = vof.copy_structure<float_t, 3>("Surface Tension Force");

        // Run interface gradient module
        tpf::modules::surface_tension<float_t> surface_tension_module;

        surface_tension_module.set_input(std::make_tuple(vof, gradients, curvature));
        surface_tension_module.set_output(surface_tension);
        surface_tension_module.set_parameters(std::make_tuple(this->Coefficient, this->Density, this->Timestep));

        surface_tension_module.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->CopyStructure(in_grid);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, surface_tension.get_name(), surface_tension.get_data(), surface_tension.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Surface Tension' failed."));

        return 0;
    }

    return 1;
}
