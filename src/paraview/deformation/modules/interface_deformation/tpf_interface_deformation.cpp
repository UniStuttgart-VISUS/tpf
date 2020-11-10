#include "tpf_interface_deformation.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/fluid_position/tpf_module_fluid_position.h"
#include "tpf_modules/interface_bending/tpf_module_interface_bending.h"
#include "tpf_modules/interface_curvature/tpf_module_interface_curvature.h"
#include "tpf_modules/interface_gradient/tpf_module_interface_gradient.h"
#include "tpf_modules/interface_stretching/tpf_module_interface_stretching.h"
#include "tpf_modules/surface_tension/tpf_module_surface_tension.h"

#include "tpf/utility/tpf_minmax.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_ghost.h"
#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <stdexcept>

vtkStandardNewMacro(tpf_interface_deformation);

tpf_interface_deformation::tpf_interface_deformation() : num_ghost_levels(0)
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_deformation::~tpf_interface_deformation() {}

int tpf_interface_deformation::FillInputPortInformation(int port, vtkInformation* info)
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

int tpf_interface_deformation::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    this->num_ghost_levels = tpf::vtk::set_ghost_levels(input_vector, output_vector, this->GetNumberOfInputPorts(),
        tpf::utility::max(this->num_ghost_levels, tpf::modules::interface_gradient<float_t>::get_num_required_ghost_levels(), tpf::modules::fluid_position<float_t>::get_num_required_ghost_levels(),
            tpf::modules::interface_curvature<float_t>::get_num_required_ghost_levels(), tpf::modules::surface_tension<float_t>::get_num_required_ghost_levels(),
            tpf::modules::interface_stretching<float_t>::get_num_required_ghost_levels(), tpf::modules::interface_bending<float_t>::get_num_required_ghost_levels()));

    return 1;
}

int tpf_interface_deformation::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto input = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, input));

        const bool has_velocity = GetInputArrayToProcess(1, input) != nullptr;

        // Create output for intermediate results
        auto gradients = vof.copy_structure<float_t, 3>("Interface Gradient");
        auto positions = vof.copy_structure<float_t, 3>("Fluid Position");
        auto curvature = vof.copy_structure<float_t, 1>("Interface Curvature");
        auto surface_tension = vof.copy_structure<float_t, 3>("Surface Tension Force");

        // Create output for temporary results
        tpf::data::polydata<float_t> position_points;

        // Create output for interface stretching
        auto stretch_area = vof.copy_structure<float_t, 1>("Stretching (area)");
        auto stretch_minimum_direction = vof.copy_structure<float_t, 3>("Stretching Direction (minimum)");
        auto stretch_maximum_direction = vof.copy_structure<float_t, 3>("Stretching Direction (maximum)");
        auto stretch_largest_direction = vof.copy_structure<float_t, 3>("Stretching Direction (largest)");

        // Create output for interface bending
        auto bend_minimum = vof.copy_structure<float_t, 1>("Bending (minimum)");
        auto bend_maximum = vof.copy_structure<float_t, 1>("Bending (maximum)");
        auto bend_largest = vof.copy_structure<float_t, 1>("Bending (absolute maximum)");
        auto bend_minimum_direction = vof.copy_structure<float_t, 3>("Bending Direction (minimum)");
        auto bend_maximum_direction = vof.copy_structure<float_t, 3>("Bending Direction (maximum)");
        auto bend_largest_direction = vof.copy_structure<float_t, 3>("Bending Direction (absolute maximum)");

        // Get parameters
        bool use_surface_tension = this->SurfaceTension != 0;

        const float_t surface_tension_coefficient = this->Coefficient;
        const float_t surface_tension_density = this->Density;

        const float_t time_step = (this->Timestep != 0.0) ? this->Timestep : tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), input);

        // Sanity check
        if (!has_velocity && !use_surface_tension)
        {
            tpf::log::warning_message(__tpf_warning_message("No velocity field specified. Using surface tension force instead."));

            use_surface_tension = true;
        }

        // Connect and run modules
        if (use_surface_tension)
        {
            // Set input, output and parameters
            tpf::modules::interface_gradient<float_t> interface_gradient_module;
            tpf::modules::fluid_position<float_t> fluid_position_module;
            tpf::modules::interface_curvature<float_t> interface_curvature_module;
            tpf::modules::surface_tension<float_t> surface_tension_module;
            tpf::modules::interface_stretching<float_t> interface_stretching_module;
            tpf::modules::interface_bending<float_t> interface_bending_module;

            interface_gradient_module.set_input(vof);
            interface_gradient_module.set_output(gradients);

            fluid_position_module.set_input(vof, gradients);
            fluid_position_module.set_parameters(tpf::modules::fluid_position_aux::position_t::INTERFACE);
            fluid_position_module.set_output(positions, position_points);

            interface_curvature_module.set_input(vof, gradients, positions);
            interface_curvature_module.set_output(curvature);

            surface_tension_module.set_input(vof, gradients, curvature);
            surface_tension_module.set_parameters(surface_tension_coefficient, surface_tension_density, time_step);
            surface_tension_module.set_output(surface_tension);

            interface_stretching_module.set_input(vof, gradients, positions, surface_tension);
            interface_stretching_module.set_parameters(time_step);
            interface_stretching_module.set_output(stretch_area, stretch_maximum_direction, stretch_minimum_direction, stretch_largest_direction);

            interface_bending_module.set_input(vof, gradients, positions, surface_tension);
            interface_bending_module.set_parameters(time_step);
            interface_bending_module.set_output(bend_minimum, bend_maximum, bend_largest, bend_minimum_direction, bend_maximum_direction, bend_largest_direction);

            // Run modules
            interface_gradient_module.run();
            fluid_position_module.run();
            interface_curvature_module.run();
            surface_tension_module.run();

            interface_stretching_module.run();
            interface_bending_module.run();
        }
        else
        {
            // Set input, output and parameters
            tpf::modules::interface_gradient<float_t> interface_gradient_module;
            tpf::modules::fluid_position<float_t> fluid_position_module;
            tpf::modules::interface_stretching<float_t> interface_stretching_module;
            tpf::modules::interface_bending<float_t> interface_bending_module;

            interface_gradient_module.set_input(vof);
            interface_gradient_module.set_output(gradients);

            fluid_position_module.set_input(vof, gradients);
            fluid_position_module.set_parameters(tpf::modules::fluid_position_aux::position_t::CELL_CENTER);
            fluid_position_module.set_output(positions, position_points);

            const auto velocity = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, input));

            interface_stretching_module.set_input(vof, gradients, positions, velocity);
            interface_stretching_module.set_parameters(time_step);
            interface_stretching_module.set_output(stretch_area, stretch_maximum_direction, stretch_minimum_direction, stretch_largest_direction);

            interface_bending_module.set_input(vof, gradients, positions, velocity);
            interface_bending_module.set_parameters(time_step);
            interface_bending_module.set_output(bend_minimum, bend_maximum, bend_largest, bend_minimum_direction, bend_maximum_direction, bend_largest_direction);

            // Run modules
            interface_gradient_module.run();
            fluid_position_module.run();

            interface_stretching_module.run();
            interface_bending_module.run();
        }

        // Set output
        auto output = vtkRectilinearGrid::GetData(output_vector);
        output->CopyStructure(input);

#ifdef __tpf_detailed
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, gradients.get_name(), gradients.get_data(), gradients.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, positions.get_name(), positions.get_data(), positions.get_num_components());

        if (use_surface_tension)
        {
            tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, curvature.get_name(), curvature.get_data(), curvature.get_num_components());
            tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, surface_tension.get_name(), surface_tension.get_data(), surface_tension.get_num_components());
        }
#endif

        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretch_area.get_name(), stretch_area.get_data(), stretch_area.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretch_minimum_direction.get_name(), stretch_minimum_direction.get_data(), stretch_minimum_direction.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretch_maximum_direction.get_name(), stretch_maximum_direction.get_data(), stretch_maximum_direction.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, stretch_largest_direction.get_name(), stretch_largest_direction.get_data(), stretch_largest_direction.get_num_components());

        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_minimum.get_name(), bend_minimum.get_data(), bend_minimum.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_maximum.get_name(), bend_maximum.get_data(), bend_maximum.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_largest.get_name(), bend_largest.get_data(), bend_largest.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_minimum_direction.get_name(), bend_minimum_direction.get_data(), bend_minimum_direction.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_maximum_direction.get_name(), bend_maximum_direction.get_data(), bend_maximum_direction.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, bend_largest_direction.get_name(), bend_largest_direction.get_data(), bend_largest_direction.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Deformation' failed."));

        return 0;
    }
    
    return 1;
}
