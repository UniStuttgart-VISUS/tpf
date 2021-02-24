#include "tpf_fs3d_reader.h"

#include "tpf/data/tpf_grid.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/fs3d_reader/tpf_module_fs3d_reader.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

vtkStandardNewMacro(tpf_fs3d_reader);

tpf_fs3d_reader::tpf_fs3d_reader()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);

    this->FileName = new char[1024];

    this->fs3d_reader = std::make_shared<tpf::modules::fs3d_reader<float_t>>();
}

tpf_fs3d_reader::~tpf_fs3d_reader()
{
    delete[] this->FileName;
}

int tpf_fs3d_reader::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    try
    {
        // Get information object
        vtkInformation* output_info = output_vector->GetInformationObject(0);

        // Read extent information
        auto& fs3d_reader = *std::static_pointer_cast<tpf::modules::fs3d_reader<float_t>>(this->fs3d_reader);

        fs3d_reader.set_parameters(this->FileName);
        auto info = fs3d_reader.provide_information();

        // Set data provided
        this->num_components = info.datasets[0].second;

        // Set time information
        this->timesteps = info.timesteps;
        double timerange[2] = { this->timesteps[0], this->timesteps[this->timesteps.size() - 1] };

        output_info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), this->timesteps.data(), static_cast<int>(this->timesteps.size()));
        output_info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timerange, 2);

        // Set extent
        int extent[6];
        extent[0] = static_cast<int>(info.extent[0].first);
        extent[1] = static_cast<int>(info.extent[0].second + 1);
        extent[2] = static_cast<int>(info.extent[1].first);
        extent[3] = static_cast<int>(info.extent[1].second + 1);
        extent[4] = static_cast<int>(info.extent[2].first);
        extent[5] = static_cast<int>(info.extent[2].second + 1);

        output_info->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

        // Can produce sub extents
        output_info->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Information request to plugin 'FS3D Reader' failed."));

        return 0;
    }

    return 1;
}

int tpf_fs3d_reader::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    try
    {
        // Get information
        vtkInformation *output_info = output_vector->GetInformationObject(0);
        auto output = vtkRectilinearGrid::GetData(output_vector);

        // Get time step
        std::size_t timestep = 0;

        if (output_info->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
        {
            const auto requested_time = output_info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
            timestep = find_timestep(requested_time);

            output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), this->timesteps[timestep]);
        }

        // Get (sub-)extent
        int subextent[6];
        output_info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), subextent);

        tpf::data::extent_t extent(3);
        extent[0] = std::make_pair(static_cast<std::size_t>(subextent[0]), static_cast<std::size_t>(subextent[1] - 1));
        extent[1] = std::make_pair(static_cast<std::size_t>(subextent[2]), static_cast<std::size_t>(subextent[3] - 1));
        extent[2] = std::make_pair(static_cast<std::size_t>(subextent[4]), static_cast<std::size_t>(subextent[5] - 1));

        // Create output data
        tpf::data::grid<float_t, float_t, 3, 1> scalar_data("Data", extent);
        tpf::data::grid<float_t, float_t, 3, 3> vector_data("Data", extent);

        using output_t = tpf::modules::fs3d_reader_aux::scalar_or_vector<float_t>;
        output_t output_data;

        // Run interface gradient module
        auto& fs3d_reader = *std::dynamic_pointer_cast<tpf::modules::fs3d_reader<float_t>>(this->fs3d_reader);

        fs3d_reader.set_run_parameters(timestep, extent);

        if (this->num_components == 1)
        {
            output_data.tag = output_t::SCALAR;
            output_data.scalar_field = &scalar_data;

            fs3d_reader.set_output(output_data);
        }
        else if (this->num_components == 3)
        {
            output_data.tag = output_t::VECTOR;
            output_data.vector_field = &vector_data;

            fs3d_reader.set_output(output_data);
        }
        else
        {
            throw std::runtime_error(__tpf_error_message("Invalid number of components from FS3D file."));
        }

        fs3d_reader.run();

        // Set output
        if (this->num_components == 1)
        {
            tpf::vtk::set_grid_information<float_t>(output, scalar_data.get_node_coordinates(), scalar_data.get_extent());
            tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, scalar_data.get_name(), scalar_data.get_data(), scalar_data.get_num_components());
        }
        else
        {
            tpf::vtk::set_grid_information<float_t>(output, vector_data.get_node_coordinates(), vector_data.get_extent());
            tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, vector_data.get_name(), vector_data.get_data(), vector_data.get_num_components());
        }
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'FS3D Reader' failed."));

        return 0;
    }

    return 1;
}

std::size_t tpf_fs3d_reader::find_timestep(const double time)
{
    double min_dist = std::numeric_limits<double>::max();
    std::size_t min = 0;

    for (std::size_t t = 0; t < this->timesteps.size(); ++t)
    {
        const double dist = std::abs(this->timesteps[t] - time);

        if (dist < min_dist)
        {
            min_dist = dist;
            min = t;
        }
    }

    return min;
}
