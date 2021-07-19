#include "tpf_correct_vof.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/correct_vof/tpf_module_correct_vof.h"

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

vtkStandardNewMacro(tpf_correct_vof);

tpf_correct_vof::tpf_correct_vof() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_correct_vof::~tpf_correct_vof() {}

int tpf_correct_vof::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_correct_vof::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::correct_vof<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_correct_vof::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_vof = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_vof, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_vof));

        // Create output data
        auto new_vof = vof.copy_structure<float_t, 1>(vof.get_name());

        // Run interface gradient module
        tpf::modules::correct_vof<float_t> correct_vof;

        correct_vof.set_input(vof);
        correct_vof.set_output(new_vof);

        correct_vof.run();

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);

        output->ShallowCopy(in_vof);
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, new_vof.get_name(), new_vof.get_data(), new_vof.get_num_components());
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Correct VOF' failed."));

        return 0;
    }

    return 1;
}
