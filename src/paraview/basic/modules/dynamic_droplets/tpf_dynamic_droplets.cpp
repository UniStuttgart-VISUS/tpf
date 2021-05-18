#include "tpf_dynamic_droplets.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/log/tpf_log.h"

#include "tpf_modules/dynamic_droplets/tpf_module_dynamic_droplets.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_log.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
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
        /// <param name="translation_name">Name of the input translation array</param>
        /// <param name="rotation_name">Name of the input rotation array</param>
        /// <param name="radius_name">Name of the input radius array</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps,
            vtkInformation* request, vtkAlgorithm* droplets_alg, std::string translation_name,
            std::string rotation_name, std::string radius_name)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->droplets_alg = droplets_alg;

            this->translation_name = translation_name;
            this->rotation_name = rotation_name;
            this->radius_name = radius_name;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, droplets]</returns>
        std::pair<float_t, tpf::data::polydata<float_t>> operator()()
        {
            // Check if next time step is available
            ++this->time_offset;

            if (this->timesteps.size() > this->time_offset)
            {
                // Request update
                request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                this->droplets_alg->Update(request);

                // Get data
                auto in_droplets = vtkPolyData::SafeDownCast(this->droplets_alg->GetOutputDataObject(0));

                auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
                    tpf::data::data_information<float_t, 3> { this->translation_name, tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 3> { this->rotation_name, tpf::data::topology_t::POINT_DATA },
                    tpf::data::data_information<float_t, 1> { this->radius_name, tpf::data::topology_t::POINT_DATA });

                // Set time step delta and return
                float_t timestep_delta = this->timesteps[this->time_offset] - this->timesteps[this->time_offset - 1];

                return std::make_pair(timestep_delta, droplets);
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

        /// Array names
        std::string translation_name;
        std::string rotation_name;
        std::string radius_name;
    };
}

vtkStandardNewMacro(tpf_dynamic_droplets);

tpf_dynamic_droplets::tpf_dynamic_droplets()
{
    this->SetNumberOfInputPorts(1);
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

        auto translation_arr = GetInputArrayToProcess(0, in_droplets);
        auto rotation_arr = GetInputArrayToProcess(1, in_droplets);
        auto radius_arr = GetInputArrayToProcess(2, in_droplets);

        if (translation_arr == nullptr || rotation_arr == nullptr || radius_arr == nullptr)
        {
            throw std::runtime_error(__tpf_error_message("Not all input arrays are provided."));
        }

        const auto droplets = tpf::vtk::get_polydata<float_t>(in_droplets,
            tpf::data::data_information<float_t, 3>{ translation_arr->GetName(), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 3>{ rotation_arr->GetName(), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<float_t, 1>{ radius_arr->GetName(), tpf::data::topology_t::POINT_DATA });

        // Create output data
        tpf::data::polydata<float_t> paths, axes, ribbons;

        if (this->ShowTranslationPaths || this->ShowRotationAxes || this->ShowRotationRibbons)
        {
            // Get time step
            const auto timestep = tpf::vtk::get_timestep_delta<float_t>(input_vector[0]->GetInformationObject(0), in_droplets);

            // Set input, output and parameters
            tpf::modules::dynamic_droplets<float_t> dynamic_droplets_module;

            dynamic_droplets_module.set_input(droplets);
            dynamic_droplets_module.set_output(
                tpf::modules::set_or_default(paths, this->ShowTranslationPaths),
                tpf::modules::set_or_default(axes, this->ShowRotationAxes),
                tpf::modules::set_or_default(ribbons, this->ShowRotationRibbons));

            dynamic_droplets_module.set_parameters(static_cast<std::size_t>(this->NumTimeSteps), static_cast<float_t>(timestep),
                static_cast<float_t>(this->RotationAxisScale), this->RotationAxisTranslation == 1,
                translation_arr->GetName(), rotation_arr->GetName(), radius_arr->GetName());

            // Create call back for loading the next time frame if requested
            const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_droplets).second;
            const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

            data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0),
                translation_arr->GetName(), rotation_arr->GetName(), radius_arr->GetName());

            dynamic_droplets_module.set_callbacks(&call_back_loader);

            // Run module
            dynamic_droplets_module.run();

            // Set output
            auto output = vtkMultiBlockDataSet::GetData(output_vector);

            auto save_lines = [](const tpf::data::polydata<float_t>& lines, vtkMultiBlockDataSet* output, const unsigned int index, const std::string& name)
            {
                auto mesh = tpf::vtk::create_mesh(lines.get_geometry());

                auto block = vtkSmartPointer<vtkPolyData>::New();
                block->SetPoints(mesh.points);
                block->SetLines(mesh.lines);

                tpf::vtk::set_data<std::size_t>(block, tpf::data::topology_t::CELL_DATA, *lines.template get_cell_data_as<std::size_t, 1>("Time ID"));

                output->SetBlock(index, block);
                output->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), name);
            };

            unsigned int block_index = 0;

            if (this->ShowTranslationPaths)
            {
                save_lines(paths, output, block_index++, "translation paths");
            }

            if (this->ShowRotationAxes)
            {
                save_lines(axes, output, block_index++, "rotation axes");
            }

            if (this->ShowRotationRibbons)
            {
                auto mesh = tpf::vtk::create_mesh(ribbons.get_geometry());

                auto block = vtkSmartPointer<vtkPolyData>::New();
                block->SetPoints(mesh.points);
                block->SetPolys(mesh.polygonals);

                tpf::vtk::set_data<std::size_t>(block, tpf::data::topology_t::POINT_DATA, *ribbons.template get_point_data_as<std::size_t, 1>("Time ID"));

                output->SetBlock(block_index, block);
                output->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "rotation ribbons");

                ++block_index;
            }

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
