#include "tpf_interface_deformation_glyph.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/interface_deformation_glyph/tpf_module_interface_deformation_glyph.h"

#include "tpf/vtk/tpf_data.h"
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

#include "vtkPolyDataNormals.h"

#include <functional>
#include <optional>
#include <stdexcept>

vtkStandardNewMacro(tpf_interface_deformation_glyph);

tpf_interface_deformation_glyph::tpf_interface_deformation_glyph()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(3);
}

tpf_interface_deformation_glyph::~tpf_interface_deformation_glyph() {}

int tpf_interface_deformation_glyph::ProcessRequest(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Create an output object of the correct type.
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
        return this->RequestDataObject(request, input_vector, output_vector);
    }

    // Generate the data
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
        return this->RequestInformation(request, input_vector, output_vector);
    }

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request, input_vector, output_vector);
    }

    if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
        return this->RequestUpdateExtent(request, input_vector, output_vector);
    }

    return this->Superclass::ProcessRequest(request, input_vector, output_vector);
}

int tpf_interface_deformation_glyph::RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    auto create_or_get_data_object = [&output_vector, this](int index)
    {
        auto output = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(index)->Get(vtkDataObject::DATA_OBJECT()));

        if (!output)
        {
            output = vtkPolyData::New();
            output_vector->GetInformationObject(index)->Set(vtkDataObject::DATA_OBJECT(), output);
            this->GetOutputPortInformation(index)->Set(vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
        }
    };

    create_or_get_data_object(0);
    create_or_get_data_object(1);
    create_or_get_data_object(2);

    return 1;
}

int tpf_interface_deformation_glyph::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }

    return 0;
}

int tpf_interface_deformation_glyph::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }
    else if (port == 1)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }
    else if (port == 2)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }

    return 1;
}

int tpf_interface_deformation_glyph::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int tpf_interface_deformation_glyph::RequestUpdateExtent(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int tpf_interface_deformation_glyph::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // Get input data
        auto input = vtkRectilinearGrid::GetData(input_vector[0]);
        const auto vof = tpf::vtk::get_grid<float_t, float_t, 3, 1>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(0, input));
        const auto positions = tpf::vtk::get_grid<float_t, float_t, 3, 3>(input, tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(1, input));

        auto get_opt_scalar = [this, &input](int index) -> std::optional<tpf::data::grid<float_t, float_t, 3, 1>>
        {
            if (GetInputArrayToProcess(index, input) != nullptr)
            {
                return std::make_optional(tpf::vtk::get_grid<float_t, float_t, 3, 1>(input,
                    tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(index, input)));
            }

            return std::nullopt;
        };

        auto get_opt_vector = [this, &input](int index) -> std::optional<tpf::data::grid<float_t, float_t, 3, 3>>
        {
            if (GetInputArrayToProcess(index, input) != nullptr)
            {
                return std::make_optional(tpf::vtk::get_grid<float_t, float_t, 3, 3>(input,
                    tpf::data::topology_t::CELL_DATA, GetInputArrayToProcess(index, input)));
            }

            return std::nullopt;
        };

        const bool velocity_glyph = output_vector->GetInformationObject(0) != nullptr && this->VelocityGlyph
            && GetInputArrayToProcess(3, input) != nullptr;
        const bool stretching_glyph = output_vector->GetInformationObject(1) != nullptr && this->StretchingGlyph
            && GetInputArrayToProcess(2, input) != nullptr && GetInputArrayToProcess(4, input) != nullptr
            && GetInputArrayToProcess(5, input) != nullptr;
        const bool bending_glyph = output_vector->GetInformationObject(2) != nullptr && this->BendingGlyph
            && GetInputArrayToProcess(2, input) != nullptr && GetInputArrayToProcess(6, input) != nullptr
            && GetInputArrayToProcess(7, input) != nullptr && GetInputArrayToProcess(8, input) != nullptr
            && GetInputArrayToProcess(9, input) != nullptr && GetInputArrayToProcess(10, input) != nullptr;

        const auto gradients = get_opt_vector(2);
        const auto velocities = get_opt_vector(3);
        const auto stretching_min = get_opt_vector(4);
        const auto stretching_max = get_opt_vector(5);
        const auto bending_min = get_opt_scalar(6);
        const auto bending_max = get_opt_scalar(7);
        const auto bending_direction_min = get_opt_vector(8);
        const auto bending_direction_max = get_opt_vector(9);
        const auto bending_polynomial = get_opt_vector(10);

        const float_t time_step = (this->Timestep != 0.0) ? this->Timestep : tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), input);

        // Create output data
        tpf::data::polydata<float_t> velocity_glyphs, stretching_glyphs, bending_glyphs;

        // Set parameters
        tpf::modules::interface_deformation_glyph_aux::velocity_params_t<float_t> velocity_parameters;
        velocity_parameters.size = static_cast<tpf::modules::interface_deformation_glyph_aux::arrow_size_t>(this->ArrowSize);
        velocity_parameters.scalar = this->ArrowScalar;
        velocity_parameters.fixed_scalar = this->ArrowFixedScalar;
        velocity_parameters.resolution = this->ArrowResolution;
        velocity_parameters.shaft_length = this->ShaftTipRatio;
        velocity_parameters.shaft_thickness = this->ArrowThickness;

        tpf::modules::interface_deformation_glyph_aux::stretching_params_t<float_t> stretching_parameters;
        stretching_parameters.exponent = this->StretchingExponent;
        stretching_parameters.size_scalar = this->StretchingSizeScalar;
        stretching_parameters.hole_radius = this->StretchingHoleRadius;
        stretching_parameters.offset = this->StretchingOffset;
        stretching_parameters.disc_resolution = this->StretchingDiscResolution;
        stretching_parameters.disc_bending = this->StretchingDiscBending;
        stretching_parameters.show_strip = this->StretchingShowStrips;
        stretching_parameters.strip_size = this->StretchingStripSize;
        stretching_parameters.show_reference = this->StretchingShowReference;
        stretching_parameters.reference_size = this->StretchingReferenceSize;
        stretching_parameters.z_offset = this->StretchingZOffset;

        tpf::modules::interface_deformation_glyph_aux::bending_params_t<float_t> bending_parameters;
        bending_parameters.scalar = this->BendingScalar;
        bending_parameters.size_scalar = this->BendingSizeScalar;
        bending_parameters.offset = this->BendingOffset;
        bending_parameters.disc_resolution = this->BendingDiscResolution;
        bending_parameters.polynomial_resolution = this->PolynomialResolution;
        bending_parameters.show_strip = this->BendingShowStrips;
        bending_parameters.strip_size = this->BendingStripSize;
        bending_parameters.z_offset = this->BendingZOffset;

        // Run interface gradient module
        tpf::modules::interface_deformation_glyph<float_t> interface_deformation_glyph;

        interface_deformation_glyph.set_input(vof, positions, gradients, velocities, stretching_min, stretching_max,
            bending_min, bending_max, bending_direction_min, bending_direction_max, bending_polynomial);
        interface_deformation_glyph.set_output(velocity_glyphs, stretching_glyphs, bending_glyphs);
        interface_deformation_glyph.set_parameters(velocity_glyph, stretching_glyph, bending_glyph, time_step,
            velocity_parameters, stretching_parameters, bending_parameters);

        interface_deformation_glyph.run();

        // Set output
        if (velocity_glyph)
        {
            auto output_velocity_glyphs = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));

            tpf::vtk::set_polydata(output_velocity_glyphs, velocity_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Velocity Magnitude", tpf::data::topology_t::OBJECT_DATA });

            auto normal_filter = vtkPolyDataNormals::New();
            normal_filter->SetInputDataObject(output_velocity_glyphs);
            normal_filter->SetFeatureAngle(180.0);
            normal_filter->Update();

            output_velocity_glyphs->ShallowCopy(normal_filter->GetOutput());
        }

        if (stretching_glyph)
        {
            auto output_stretching_glyphs = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));

            tpf::vtk::set_polydata(output_stretching_glyphs, stretching_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Stretching", tpf::data::topology_t::OBJECT_DATA });

            auto normal_filter = vtkPolyDataNormals::New();
            normal_filter->SetInputDataObject(output_stretching_glyphs);
            normal_filter->SetFeatureAngle(180.0);
            normal_filter->Update();

            output_stretching_glyphs->ShallowCopy(normal_filter->GetOutput());
        }

        if (bending_glyph)
        {
            auto output_bending_glyphs = vtkPolyData::SafeDownCast(output_vector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));

            tpf::vtk::set_polydata(output_bending_glyphs, bending_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Bending", tpf::data::topology_t::OBJECT_DATA });

            auto normal_filter = vtkPolyDataNormals::New();
            normal_filter->SetInputDataObject(output_bending_glyphs);
            normal_filter->SetFeatureAngle(180.0);
            normal_filter->Update();

            output_bending_glyphs->ShallowCopy(normal_filter->GetOutput());
        }
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Deformation Glyph' failed."));

        return 0;
    }

    return 1;
}
