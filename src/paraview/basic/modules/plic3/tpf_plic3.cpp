#include "tpf_plic3.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/plic3/tpf_module_plic3.h"

#include "tpf/vtk/tpf_grid.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkCellData.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"

#include <stdexcept>
#include <tuple>

vtkStandardNewMacro(tpf_plic3);

tpf_plic3::tpf_plic3()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_plic3::~tpf_plic3()
{
}

int tpf_plic3::FillInputPortInformation(int port, vtkInformation *info)
{
    if (!this->Superclass::FillInputPortInformation(port, info)) {
        return 0;
    }
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
}

int tpf_plic3::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    try
    {
        // get the grid
        vtkRectilinearGrid *grid = vtkRectilinearGrid::GetData(input_vector[0]);

        // data arrays from GUI selection
        vtkDataArray *data_array_f = this->GetInputArrayToProcess(0, input_vector);
        vtkDataArray *data_array_f3 = this->GetInputArrayToProcess(1, input_vector);
        vtkDataArray *data_array_f_norm_3ph = this->GetInputArrayToProcess(2, input_vector);

        // create TPF grids
        // TODO selected data_arrays are now only used to get names, perhaps it is more efficient to use the data_arrays directly with a new get_grid function.
        const auto f = tpf::vtk::get_grid<float_t, float_t, 3, 1>(grid, tpf::data::topology_t::CELL_DATA, data_array_f->GetName());
        const auto f3 = tpf::vtk::get_grid<float_t, float_t, 3, 1>(grid, tpf::data::topology_t::CELL_DATA, data_array_f3->GetName());
        const auto f_norm_3ph = tpf::vtk::get_grid<float_t, float_t, 3, 3>(grid, tpf::data::topology_t::CELL_DATA, data_array_f_norm_3ph->GetName());

        // Create output data
        tpf::data::polydata<float_t> plic3_interface;

        // Run PLIC module
        tpf::modules::plic3<float_t> plic3_module;
        plic3_module.set_input(f, f3, f_norm_3ph);
        plic3_module.set_output(plic3_interface);
        plic3_module.set_parameters(this->Epsilon, this->ErrorMarginPlic, this->ErrorMarginPlic3, this->NumIterationsPlic,
            this->NumIterationsPlic3, this->Perturbation, static_cast<bool>(this->HideErrorCubes));
        plic3_module.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        tpf::vtk::set_polydata(output, plic3_interface,
                tpf::data::data_information<int, 3>{ "coords", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<int, 1>{ "type", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<float_t, 3>{ "grad", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<float_t, 3>{ "grad3", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<float_t, 1>{ "f", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<float_t, 1>{ "f3", tpf::data::topology_t::OBJECT_DATA },
                tpf::data::data_information<float_t, 3>{ "fnorm3ph", tpf::data::topology_t::OBJECT_DATA });
    }
    catch (const std::runtime_error& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'PLIC3' failed."));

        return 0;
    }

    return 1;
}
