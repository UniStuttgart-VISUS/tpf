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

#include <functional>
#include <optional>
#include <stdexcept>

vtkStandardNewMacro(tpf_interface_deformation_glyph);

tpf_interface_deformation_glyph::tpf_interface_deformation_glyph()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_interface_deformation_glyph::~tpf_interface_deformation_glyph() {}

int tpf_interface_deformation_glyph::FillInputPortInformation(int port, vtkInformation* info)
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

        const bool velocity_glyph = this->VelocityGlyph && GetInputArrayToProcess(2, input) != nullptr;
        const bool stretching_glyph = this->StretchingGlyph && GetInputArrayToProcess(3, input) != nullptr && GetInputArrayToProcess(4, input) != nullptr;
        const bool bending_glyph = this->BendingGlyph && GetInputArrayToProcess(5, input) != nullptr && GetInputArrayToProcess(6, input) != nullptr
            && GetInputArrayToProcess(7, input) != nullptr && GetInputArrayToProcess(8, input) != nullptr;

        const auto velocities = get_opt_vector(2);
        const auto stretching_min = get_opt_vector(3);
        const auto stretching_max = get_opt_vector(4);
        const auto bending_min = get_opt_scalar(5);
        const auto bending_max = get_opt_scalar(6);
        const auto bending_direction_min = get_opt_vector(7);
        const auto bending_direction_max = get_opt_vector(8);

        const float_t time_step = (this->Timestep != 0.0) ? this->Timestep : tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), input);

        // Create output data
        tpf::data::polydata<float_t> velocity_glyphs, stretching_glyphs, bending_glyphs;

        // Run interface gradient module
        tpf::modules::interface_deformation_glyph<float_t> interface_deformation_glyph;

        interface_deformation_glyph.set_input(vof, positions, velocities, stretching_min, stretching_max,
            bending_min, bending_max, bending_direction_min, bending_direction_max);
        interface_deformation_glyph.set_output(velocity_glyphs, stretching_glyphs, bending_glyphs);
        interface_deformation_glyph.set_parameters(velocity_glyph, stretching_glyph, bending_glyph, time_step,
            static_cast<tpf::modules::interface_deformation_glyph_aux::arrow_size_t>(this->ArrowSize),
            this->ArrowScalar, this->ArrowFixedScalar, this->ArrowResolution, this->ShaftTipRatio, this->ArrowThickness,
            this->BendingDiscResolution, this->PolygonalResolution, this->BendingStripSize, this->BendingSizeScalar, this->BendingScalar);

        interface_deformation_glyph.run();

        // Set output
        auto output = vtkMultiBlockDataSet::GetData(output_vector);

        unsigned int index = 0;

        if (velocity_glyph)
        {
            auto output_velocity_glyphs = vtkPolyData::New();

            tpf::vtk::set_polydata(output_velocity_glyphs, velocity_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Velocity Magnitude", tpf::data::topology_t::OBJECT_DATA });

            output->SetBlock(index, output_velocity_glyphs);
            output->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "velocity");

            ++index;
        }

        if (stretching_glyph)
        {
            auto output_stretching_glyphs = vtkPolyData::New();

            tpf::vtk::set_polydata(output_stretching_glyphs, stretching_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Stretching", tpf::data::topology_t::OBJECT_DATA });

            output->SetBlock(index, output_stretching_glyphs);
            output->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "stretching");

            ++index;
        }

        if (bending_glyph)
        {
            auto output_bending_glyphs = vtkPolyData::New();

            tpf::vtk::set_polydata(output_bending_glyphs, bending_glyphs,
                tpf::data::data_information<float_t, 1, 1>{ "Bending", tpf::data::topology_t::OBJECT_DATA });

            output->SetBlock(index, output_bending_glyphs);
            output->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), "bending");

            ++index;
        }
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Interface Deformation Glyph' failed."));

        return 0;
    }

    return 1;
}
