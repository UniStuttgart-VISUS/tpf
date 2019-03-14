#include "tpf_ghost.h"

#include "../log/tpf_log.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <stdexcept>

namespace tpf
{
    namespace vtk
    {
        inline std::size_t set_ghost_levels(vtkInformationVector **input_vector, vtkInformationVector *output_vector,
            const std::size_t num_input_vectors, const std::size_t num_required_ghost_levels)
        {
            try
            {
                // Get number of ghost cells required by filters downstream
                const std::size_t num_ghost_cells = std::max(num_required_ghost_levels, static_cast<std::size_t>(
                    output_vector->GetInformationObject(0)->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS())));

                // Set number of ghost cells required
                for (std::size_t i = 0; i < num_input_vectors; ++i)
                {
                    if (input_vector[i]->GetInformationObject(0) != nullptr)
                    {
                        input_vector[i]->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), static_cast<int>(num_ghost_cells));
                    }
                }

                return num_ghost_cells;
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Failed to set ghost levels."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Failed to set ghost levels."));
            }
        }
    }
}