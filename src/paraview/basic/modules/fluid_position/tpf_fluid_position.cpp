#include "tpf_fluid_position.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/fluid_position/tpf_module_fluid_position.h"

#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <stdexcept>

vtkStandardNewMacro(tpf_fluid_position);

tpf_fluid_position::tpf_fluid_position() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_fluid_position::~tpf_fluid_position() {}

int tpf_fluid_position::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_fluid_position::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::fluid_position<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_fluid_position::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));

        const bool has_gradients = GetInputArrayToProcess(1, in_grid) != nullptr;

        // Create output data
        auto positions_grid = vof.copy_structure<float_t, 3>("Fluid Positions");

        tpf::data::polydata<float_t> positions_points;

        // Run interface gradient module
        if (!has_gradients)
        {
            tpf::modules::fluid_position<float_t> fluid_position;

            fluid_position.set_input(vof, std::nullopt);
            fluid_position.set_output(positions_grid, positions_points);
            fluid_position.set_parameters(static_cast<tpf::modules::fluid_position_aux::position_t>(this->PositionType));

            fluid_position.run();
        }
        else
        {
            const auto gradients = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, in_grid));

            tpf::modules::fluid_position<float_t> fluid_position;

            fluid_position.set_input(vof, gradients);
            fluid_position.set_output(positions_grid, positions_points);
            fluid_position.set_parameters(static_cast<tpf::modules::fluid_position_aux::position_t>(this->PositionType));

            fluid_position.run();
        }
        
        // Set output grid
        auto output = vtkMultiBlockDataSet::GetData(output_vector);
        auto output_grid = vtkSmartPointer<vtkRectilinearGrid>::New();

        output_grid->ShallowCopy(in_grid);

        tpf::vtk::set_data<float_t>(output_grid, tpf::data::topology_t::CELL_DATA, positions_grid.get_name(), positions_grid.get_data(), positions_grid.get_num_components());

        output->SetBlock(0u, output_grid);
        output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), "grid");

        // Set output points
        auto output_positions = vtkSmartPointer<vtkPolyData>::New();
        tpf::vtk::set_polydata(output_positions, positions_points);

        output->SetBlock(1u, output_positions);
        output->GetMetaData(1u)->Set(vtkCompositeDataSet::NAME(), "positions");
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Fluid Position' failed."));

        return 0;
    }

    return 1;
}
