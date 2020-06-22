#include "tpf_fs3d_writer.h"

#include "tpf/data/tpf_data_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/fs3d_writer/tpf_module_fs3d_writer.h"

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
#include <string>

vtkStandardNewMacro(tpf_fs3d_writer);

tpf_fs3d_writer::tpf_fs3d_writer() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    this->FileName = new char[1024];
    this->Name = new char[1024];
    this->Unit = new char[1024];
    this->GridUnit = new char[1024];
}

tpf_fs3d_writer::~tpf_fs3d_writer()
{
    delete[] this->FileName;
    delete[] this->Name;
    delete[] this->Unit;
    delete[] this->GridUnit;
}

int tpf_fs3d_writer::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_fs3d_writer::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::fs3d_writer<float_t, 3>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_fs3d_writer::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);

        const auto timestep = tpf::vtk::get_timestep<double>(input_vector[0]->GetInformationObject(0), in_grid);

        // Write VOF field if present
        bool found = false;

        try
        {
            const auto vof_grid = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));
            found = true;

            tpf::modules::fs3d_writer<float_t, 1> fs3d_writer;

            fs3d_writer.set_input(vof_grid);
            fs3d_writer.set_parameters(std::string(this->FileName));
            fs3d_writer.set_run_parameters(static_cast<std::size_t>(timestep.second),
                vof_grid.get_extent(), timestep.first, std::string(this->Name), std::string(this->Unit), std::string(this->GridUnit));

            fs3d_writer.run();
        }
        catch (const std::runtime_error& ex)
        {
            if (found)
            {
                throw ex;
            }
        }

        // Write velocity field if present
        found = false;
        
        try
        {
            const auto velocity_grid = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));
            found = true;

            tpf::modules::fs3d_writer<float_t, 3> fs3d_writer;

            fs3d_writer.set_input(velocity_grid);
            fs3d_writer.set_parameters(std::string(this->FileName));
            fs3d_writer.set_run_parameters(static_cast<std::size_t>(timestep.second),
                velocity_grid.get_extent(), timestep.first, std::string(this->Name), std::string(this->Unit), std::string(this->GridUnit));

            fs3d_writer.run();
        }
        catch (const std::runtime_error& ex)
        {
            if (found)
            {
                throw ex;
            }
        }

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);
        output->ShallowCopy(in_grid);
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Correct VOF' failed."));

        return 0;
    }

    return 1;
}
