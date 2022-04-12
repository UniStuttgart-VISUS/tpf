#include "tpf_dynamic_droplets.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/dynamic_droplets/tpf_module_dynamic_droplets.h"

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

#include <algorithm>
#include <exception>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace
{
    /// Call back class for next time step
    /// <template name="float_t">Floating point type</template>
    template <typename float_t>
    class data_handler : public tpf::modules::dynamic_droplets_aux::request_frame_call_back<float_t>
    {
    public:
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="current_timestep">Current time step this plugin was executed in</param>
        /// <param name="timesteps">Available time steps</param>
        /// <param name="request">Original request that will be modified for requesting different time frames</param>
        /// <param name="droplets_alg">Algorithm producing droplet data</param>
        /// <param name="droplet_ids_alg">Algorithm producing droplet ID data</param>
        /// <param name="translation_name">Name of the input translation array</param>
        /// <param name="rotation_name">Name of the input rotation array</param>
        /// <param name="radius_name">Name of the input radius array</param>
        /// <param name="droplet_ids_name">Name of the input droplet ID array</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps,
            vtkInformation* request, vtkAlgorithm* droplets_alg, vtkAlgorithm* droplet_ids_alg, std::string translation_name,
            std::string rotation_name, std::string radius_name, std::string droplet_ids_name)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->droplets_alg = droplets_alg;
            this->droplet_ids_alg = droplet_ids_alg;

            this->translation_name = translation_name;
            this->rotation_name = rotation_name;
            this->radius_name = radius_name;
            this->droplet_ids_name = droplet_ids_name;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, droplets]</returns>
        std::tuple<float_t, tpf::data::polydata<float_t>, tpf::data::grid<long long, float_t, 3, 1>> operator()()
        {
            // Check if next time step is available
            ++this->time_offset;

            if (this->timesteps.size() > this->time_offset)
            {
                // Request update
                request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                this->droplets_alg->Update(request);
                this->droplet_ids_alg->Update(request);

                // Get data
                auto in_droplets = vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0));
                auto in_droplet_ids = vtkRectilinearGrid::SafeDownCast(this->droplet_ids_alg->GetOutputDataObject(0));

                auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
                    tpf::data::data_information<float_t, 3> { this->translation_name, tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 3> { this->rotation_name, tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 1> { this->radius_name, tpf::data::topology_t::POINT_DATA });

                const auto droplet_ids = tpf::vtk::get_grid<long long, float_t, 3, 1>(in_droplet_ids,
                    tpf::data::topology_t::CELL_DATA, this->droplet_ids_name);

                // Set time step delta and return
                float_t timestep_delta = this->timesteps[this->time_offset] - this->timesteps[this->time_offset - 1];

                return std::make_tuple(timestep_delta, droplets, droplet_ids);
            }

            throw std::exception();
        }

        /// <summary>
        /// Reset input algorithms to their beginning state
        /// </summary>
        void reset()
        {
            // Request update to original time step
            this->time_offset = this->original_time;

            request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

            this->droplets_alg->Update(request);
            this->droplet_ids_alg->Update(request);
        }

    private:
        /// Original timestep
        std::size_t original_time;

        /// Last requested timestep
        std::size_t time_offset;

        /// Available timesteps
        std::vector<float_t> timesteps;

        /// Request
        vtkInformation* request;

        /// Droplet algorithm
        vtkAlgorithm* droplets_alg;
        vtkAlgorithm* droplet_ids_alg;

        /// Array names
        std::string translation_name;
        std::string rotation_name;
        std::string radius_name;
        std::string droplet_ids_name;
    };

    enum class geometry_t
    {
        point,
        line,
        polygon
    };

    template <typename data_t, std::size_t components>
    void save_block(vtkMultiBlockDataSet* output, const unsigned int block_index, const tpf::data::polydata<float_t>& data,
        const geometry_t geometry, const tpf::data::topology_t topology, const std::string& data_name, const std::string& block_name)
    {
        auto block = vtkSmartPointer<vtkPolyData>::New();

        auto mesh = tpf::vtk::create_mesh(data.get_geometry());

        block->SetPoints(mesh.points);

        if (geometry == geometry_t::point)
        {
            block->SetVerts(mesh.vertices);
        }
        else if (geometry == geometry_t::line)
        {
            block->SetLines(mesh.lines);
        }
        else
        {
            block->SetPolys(mesh.polygonals);
        }

        if (topology == tpf::data::topology_t::POINT_DATA)
        {
            tpf::vtk::set_data<data_t>(block, topology, *data.template get_point_data_as<data_t, components>(data_name));
        }
        else if (topology == tpf::data::topology_t::CELL_DATA)
        {
            tpf::vtk::set_data<data_t>(block, topology, *data.template get_cell_data_as<data_t, components>(data_name));
        }

        output->SetBlock(block_index, block);
        output->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), block_name);
    }
}

vtkStandardNewMacro(tpf_dynamic_droplets);

tpf_dynamic_droplets::tpf_dynamic_droplets()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

tpf_dynamic_droplets::~tpf_dynamic_droplets() {}

int tpf_dynamic_droplets::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int tpf_dynamic_droplets::RequestUpdateExtent(vtkInformation *vtkNotUsed(request), vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    return 1;
}

int tpf_dynamic_droplets::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data
        auto in_droplets = vtkPolyData::GetData(input_vector[0]);
        auto in_droplet_ids = vtkRectilinearGrid::GetData(input_vector[1]);

        auto translation_arr = GetInputArrayToProcess(0, in_droplets);
        auto rotation_arr = GetInputArrayToProcess(1, in_droplets);
        auto radius_arr = GetInputArrayToProcess(2, in_droplets);
        auto ids_arr = GetInputArrayToProcess(3, in_droplet_ids);

        if (translation_arr == nullptr || rotation_arr == nullptr || radius_arr == nullptr || ids_arr == nullptr)
        {
            throw std::runtime_error(__tpf_error_message("Not all input arrays are provided."));
        }

        const auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
            tpf::data::data_information<float_t, 3>{ translation_arr->GetName(), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 3>{ rotation_arr->GetName(), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ radius_arr->GetName(), tpf::data::topology_t::POINT_DATA });

        const auto droplet_ids = tpf::vtk::get_grid<long long, float_t, 3, 1>(in_droplet_ids, tpf::data::topology_t::CELL_DATA, ids_arr);

        // Create output data
        tpf::data::polydata<float_t> droplet_tracks, summary, paths, axes, ribbons, rotation_paths, coordinate_axes;

        if (this->Compute)
        {
            // Get time step
            const auto timestep = tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), in_droplets);

            // Set input, output and parameters
            tpf::modules::dynamic_droplets<float_t> dynamic_droplets_module;

            dynamic_droplets_module.set_input(droplets, droplet_ids);
            dynamic_droplets_module.set_output(droplet_tracks, summary, paths, axes, ribbons, rotation_paths, coordinate_axes);

            dynamic_droplets_module.set_parameters(static_cast<std::size_t>(this->NumTimeSteps), static_cast<float_t>(timestep),
                this->StaticFrameOfReference == 1, static_cast<float_t>(this->RibbonSize), static_cast<float_t>(this->RibbonThickness),
                this->FixRotationAxisSize == 1, static_cast<float_t>(this->RotationAxisScale), this->RotationAxisTranslation == 1,
                translation_arr->GetName(), rotation_arr->GetName(), radius_arr->GetName());

            // Create call back for loading the next time frame if requested
            const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_droplets).second;
            const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

            data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(1, 0),
                translation_arr->GetName(), rotation_arr->GetName(), radius_arr->GetName(), ids_arr->GetName());

            dynamic_droplets_module.set_callbacks(&call_back_loader);

            // Run module
            dynamic_droplets_module.run();

            // Set output
            auto output = vtkMultiBlockDataSet::GetData(output_vector);

            unsigned int index = 0;

            save_block<std::size_t, 1>(output, index++, droplet_tracks, geometry_t::line, tpf::data::topology_t::CELL_DATA, "Topology information", "droplet tracks");
            save_block<float_t, 3>(output, index++, summary, geometry_t::point, tpf::data::topology_t::POINT_DATA, "Displacement", "translation summary");
            save_block<std::size_t, 1>(output, index++, paths, geometry_t::line, tpf::data::topology_t::POINT_DATA, "Time ID", "translation paths");
            save_block<std::size_t, 1>(output, index++, axes, geometry_t::line, tpf::data::topology_t::CELL_DATA, "Time ID", "rotation axes");
            save_block<std::size_t, 1>(output, index++, ribbons, geometry_t::polygon, tpf::data::topology_t::POINT_DATA, "Time ID", "rotation ribbons");
            save_block<std::size_t, 1>(output, index++, rotation_paths, geometry_t::line, tpf::data::topology_t::POINT_DATA, "Time ID", "rotation paths");
            save_block<float_t, 1>(output, index++, coordinate_axes, geometry_t::polygon, tpf::data::topology_t::POINT_DATA, "Time information", "coordinate axes");

            // Copy time information
            for (unsigned int i = 0; i < output->GetNumberOfBlocks(); ++i)
            {
                tpf::vtk::copy_time_information(vtkPolyData::SafeDownCast(
                    GetInputAlgorithm(0, 0)->GetOutputDataObject(0)), vtkDataSet::SafeDownCast(output->GetBlock(i)));
            }
        }
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Dynamic Droplets' failed."));

        return 0;
    }

    return 1;
}
