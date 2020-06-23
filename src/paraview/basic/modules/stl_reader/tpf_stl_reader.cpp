#include "tpf_stl_reader.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/stl_reader/tpf_module_stl_reader.h"

#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <memory>
#include <stdexcept>

vtkStandardNewMacro(tpf_stl_reader);

tpf_stl_reader::tpf_stl_reader()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);

    this->FileName = new char[1024];

    this->stl_reader = std::make_shared<tpf::modules::stl_reader<float_t>>();
}

tpf_stl_reader::~tpf_stl_reader()
{
    delete[] this->FileName;
}

int tpf_stl_reader::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    try
    {
        // Get information object
        vtkInformation* output_info = output_vector->GetInformationObject(0);

        // Read extent information
        tpf::modules::stl_reader<float_t>& stl_reader_module = *std::static_pointer_cast<tpf::modules::stl_reader<float_t>>(this->stl_reader);

        stl_reader_module.set_parameters(this->FileName);
        auto info = stl_reader_module.provide_information();

        // Set data provided
        this->num_triangles = info.datasets[0].second;
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Information request to plugin 'STL Reader' failed."));

        return 0;
    }

    return 1;
}

int tpf_stl_reader::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    try
    {
        // Get information
        vtkInformation *output_info = output_vector->GetInformationObject(0);
        auto output = vtkPolyData::GetData(output_vector);

        // Create output data
        tpf::data::polydata<float_t> triangles;

        using output_t = tpf::modules::stl_reader_aux::polydata_or_buffer<float_t>;
        output_t output_triangles{ output_t::POLYDATA, &triangles };

        // Run interface gradient module
        tpf::modules::stl_reader<float_t>& stl_reader_module = *std::static_pointer_cast<tpf::modules::stl_reader<float_t>>(this->stl_reader);

        stl_reader_module.set_additional_run_parameter(this->CalculateNormals != 0);
        stl_reader_module.set_output(output_triangles);

        stl_reader_module.run();

        // Set output
        tpf::vtk::set_polydata(output, triangles
            , tpf::data::data_information<float_t, 3>{ "Surface Normal", tpf::data::topology_t::CELL_DATA }
#if defined(__tpf_sanity_checks) && defined(__tpf_detailed)
            , tpf::data::data_information<unsigned char, 1>{ "Valid", tpf::data::topology_t::CELL_DATA }
#endif
        );
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'STL Reader' failed."));

        return 0;
    }

    return 1;
}
