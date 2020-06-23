#include "tpf_test_data.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/test_data/tpf_module_test_data.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_vtk_traits.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <stdexcept>

vtkStandardNewMacro(tpf_test_data);

tpf_test_data::tpf_test_data()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

tpf_test_data::~tpf_test_data() {}

int tpf_test_data::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    vtkInformation* output_info = output_vector->GetInformationObject(0);

    // Set extent
    int extent[6];
    extent[0] = static_cast<int>(0);
    extent[1] = static_cast<int>(this->NumCells);
    extent[2] = static_cast<int>(0);
    extent[3] = static_cast<int>(this->NumCells);
    extent[4] = static_cast<int>(0);
    extent[5] = static_cast<int>(this->NumCells);

    output_info->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

    // Cannot produce sub extents
    output_info->Set(CAN_PRODUCE_SUB_EXTENT(), 0);

    return 1;
}

int tpf_test_data::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    try
    {
        // Create output
        tpf::data::grid<float_t, float_t, 3, 1> vof;
        tpf::data::grid<float_t, float_t, 3, 3> velocities;

        // Run test data module
        tpf::modules::test_data<float_t> test_data_module;

        test_data_module.set_output(vof, velocities);
        test_data_module.set_parameters(static_cast<tpf::modules::test_data_aux::geometry_t>(this->Geometry), this->NumCells,
            this->Spinning == 1, this->Moving == 1, this->Expanding == 1, this->Rotating == 1,
            this->SpinningMagnitude, this->MovingMagnitude, this->ExpandingMagnitude, this->RotatingMagnitude);

        test_data_module.run();

        // Set output structure
        auto output = vtkRectilinearGrid::GetData(output_vector);
        
        int extent[6] = {
            static_cast<int>(vof.get_extent()[0].first),
            static_cast<int>(vof.get_extent()[0].second + 1),
            static_cast<int>(vof.get_extent()[1].first),
            static_cast<int>(vof.get_extent()[1].second + 1),
            static_cast<int>(vof.get_extent()[2].first),
            static_cast<int>(vof.get_extent()[2].second + 1) };

        output->SetExtent(extent);

        auto x_coordinates = vtkSmartPointer<typename tpf::vtk::vtk_array<float_t>::type>::New();
        auto y_coordinates = vtkSmartPointer<typename tpf::vtk::vtk_array<float_t>::type>::New();
        auto z_coordinates = vtkSmartPointer<typename tpf::vtk::vtk_array<float_t>::type>::New();

        x_coordinates->SetNumberOfComponents(1);
        y_coordinates->SetNumberOfComponents(1);
        z_coordinates->SetNumberOfComponents(1);
        
        x_coordinates->SetNumberOfValues(vof.get_node_coordinates()[0].size());
        y_coordinates->SetNumberOfValues(vof.get_node_coordinates()[1].size());
        z_coordinates->SetNumberOfValues(vof.get_node_coordinates()[2].size());

        for (std::size_t i = 0; i < static_cast<std::size_t>(x_coordinates->GetNumberOfValues()); ++i)
        {
            x_coordinates->SetValue(i, vof.get_node_coordinates()[0][i]);
        }
        for (std::size_t i = 0; i < static_cast<std::size_t>(y_coordinates->GetNumberOfValues()); ++i)
        {
            y_coordinates->SetValue(i, vof.get_node_coordinates()[1][i]);
        }
        for (std::size_t i = 0; i < static_cast<std::size_t>(z_coordinates->GetNumberOfValues()); ++i)
        {
            z_coordinates->SetValue(i, vof.get_node_coordinates()[2][i]);
        }

        output->SetXCoordinates(x_coordinates);
        output->SetYCoordinates(y_coordinates);
        output->SetZCoordinates(z_coordinates);

        // Set output data
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, "Volume of Fluid", vof.get_data(), vof.get_num_components());
        tpf::vtk::set_data<float_t>(output, tpf::data::topology_t::CELL_DATA, "Velocity", velocities.get_data(), velocities.get_num_components());
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Test Data' failed."));

        return 0;
    }

    return 1;
}
