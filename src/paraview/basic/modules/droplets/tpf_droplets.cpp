#include "tpf_droplets.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/droplets/tpf_module_droplets.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkIdTypeArray.h"

#include <algorithm>
#include <exception>
#include <functional>
#include <optional>
#include <tuple>

vtkStandardNewMacro(tpf_droplets);

tpf_droplets::tpf_droplets() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
}

tpf_droplets::~tpf_droplets() { }

int tpf_droplets::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }

    return 0;
}

int tpf_droplets::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkRectilinearGrid");
        return 1;
    }
    else if (port == 1)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }

    return 1;
}

int tpf_droplets::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        std::max(this->num_ghost_levels, tpf::modules::droplets<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_droplets::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto in_grid = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, in_grid));
        const auto positions = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, in_grid));

        const bool has_velocity = GetInputArrayToProcess(2, in_grid) != nullptr;

        // Create output data
        auto droplet_ids = vof.copy_structure<long long, 1>("Droplet ID");
        auto droplet_volumes = vof.copy_structure<float_t, 1>("Volume");

        tpf::data::polydata<float_t> droplets;

        auto output_grid = vtkRectilinearGrid::SafeDownCast(output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
        output_grid->ShallowCopy(in_grid);

        // Run droplet module
        tpf::modules::droplets<float_t> droplets_module;

        droplets_module.set_parameters(this->CalculateTranslation != 0, this->CalculateRotation != 0, this->CalculateEnergy != 0,
            this->CalculateInertia != 0, static_cast<tpf::modules::droplets_aux::scale_method_t>(this->Scale),
            static_cast<tpf::modules::droplets_aux::rotation_method_t>(this->RotationMethod));

        if (!has_velocity)
        {
            droplets_module.set_input(vof, positions, std::nullopt);
            droplets_module.set_output(droplet_ids, droplet_volumes, std::nullopt, droplets);

            droplets_module.run();
        }
        else
        {
            const auto velocities = tpf::vtk::get_grid<float_t, float_t, 3, 3>(in_grid, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(2, in_grid));

            auto global_velocities = vof.copy_structure<float_t, 3>("Global Velocity");

            droplets_module.set_input(vof, positions, velocities);
            droplets_module.set_output(droplet_ids, droplet_volumes, global_velocities, droplets);

            droplets_module.run();

            if (this->CalculateTranslation != 0 || this->CalculateRotation != 0)
            {
                tpf::vtk::set_data<float_t>(output_grid, tpf::data::topology_t::CELL_DATA, global_velocities.get_name(), global_velocities.get_data(), global_velocities.get_num_components());
            }
        }

        // Set output grid
        tpf::vtk::set_data<long long>(output_grid, tpf::data::topology_t::CELL_DATA, droplet_ids.get_name(), droplet_ids.get_data(), droplet_ids.get_num_components());
        tpf::vtk::set_data<float_t>(output_grid, tpf::data::topology_t::CELL_DATA, droplet_volumes.get_name(), droplet_volumes.get_data(), droplet_volumes.get_num_components());

        // Set output points
        auto output_positions = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
        tpf::vtk::set_polydata(output_positions, droplets, tpf::data::data_information<long long, 1>{ "Droplet ID", tpf::data::topology_t::POINT_DATA });

        if (has_velocity)
        {
            if (this->CalculateTranslation != 0)
            {
                tpf::vtk::append_polydata(output_positions, droplets,
                    tpf::data::data_information<float_t, 3>{ "Translation", tpf::data::topology_t::POINT_DATA });

                if (this->CalculateEnergy != 0)
                {
                    tpf::vtk::append_polydata(output_positions, droplets,
                        tpf::data::data_information<float_t, 1>{ "Translation Energy", tpf::data::topology_t::POINT_DATA });
                }
            }

            if (this->CalculateRotation != 0)
            {
                tpf::vtk::append_polydata(output_positions, droplets,
                    tpf::data::data_information<float_t, 3>{ "Rotation", tpf::data::topology_t::POINT_DATA });

                if (this->CalculateEnergy != 0)
                {
                    tpf::vtk::append_polydata(output_positions, droplets,
                        tpf::data::data_information<float_t, 1>{ "Rotation Energy", tpf::data::topology_t::POINT_DATA },
                        tpf::data::data_information<float_t, 1>{ "Local Energy", tpf::data::topology_t::POINT_DATA });
                }
            }

            if (this->CalculateInertia != 0)
            {
                tpf::vtk::append_polydata(output_positions, droplets,
                    tpf::data::data_information<float_t, 3, 3>{ "Inertia", tpf::data::topology_t::POINT_DATA });

#ifdef __tpf_detailed
                tpf::vtk::append_polydata(output_positions, droplets,
                    tpf::data::data_information<float_t, 3>{ "Inertia #1", tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 3>{ "Inertia #2", tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 3>{ "Inertia #3", tpf::data::topology_t::POINT_DATA });
#endif
            }

            if (this->CalculateEnergy != 0)
            {
                tpf::vtk::append_polydata(output_positions, droplets,
                    tpf::data::data_information<float_t, 1>{ "Energy", tpf::data::topology_t::POINT_DATA });
            }
        }

        tpf::vtk::append_polydata(output_positions, droplets,
            tpf::data::data_information<float_t, 1>{ "Volume", tpf::data::topology_t::POINT_DATA });

        tpf::vtk::append_polydata(output_positions, droplets,
            tpf::data::data_information<float_t, 1>{ "Radius", tpf::data::topology_t::POINT_DATA });
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Droplets' failed."));

        return 0;
    }

    return 1;
}
