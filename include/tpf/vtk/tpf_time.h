#pragma once

#include "vtkDataSet.h"
#include "vtkInformation.h"

#include <utility>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        /// <summary>
        /// Get the current time step
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <param name="info">VTK information object</param>
        /// <return>[time step, time step index]</return>
        template <typename float_t>
        std::vector<float_t> get_timesteps(vtkInformation* info);

        /// <summary>
        /// Get the current time step
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <param name="info">VTK information object</param>
        /// <param name="alg">VTK algorithm object</param>
        /// <return>[time step, time step index]</return>
        template <typename float_t>
        std::pair<float_t, std::size_t> get_timestep(vtkInformation* info, vtkDataSet* alg);

        /// <summary>
        /// Get the time step at given index
        /// </summary>
        /// <template name="T">Scalar type</template>
        /// <param name="info">VTK information object</param>
        /// <param name="timestep_index">Time step index</param>
        /// <return>Time step</return>
        template <typename float_t>
        float_t get_timestep(vtkInformation* info, std::size_t timestep_index);

        /// <summary>
        /// Is the time step at given index the first time step?
        /// </summary>
        /// <param name="info">VTK information object</param>
        /// <param name="timestep_index">Time step index</param>
        /// <return>True if it is the first time step</return>
        bool is_first_timestep(vtkInformation* info, std::size_t timestep_index);

        /// <summary>
        /// Is the time step at given index the last time step?
        /// </summary>
        /// <param name="info">VTK information object</param>
        /// <param name="timestep_index">Time step index</param>
        /// <return>True if it is the last time step</return>
        bool is_last_timestep(vtkInformation* info, std::size_t timestep_index);

        /// <summary>
        /// Get the current time step delta
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <param name="info">VTK information object</param>
        /// <param name="alg">VTK algorithm object</param>
        /// <return>Current time step delta</return>
        template <typename float_t>
        float_t get_timestep_delta(vtkInformation* info, vtkDataSet* alg);

        /// <summary>
        /// Copy time range and time values
        /// </summary>
        /// <param name="source_info">Information object, which holds the time range and time values</param>
        /// <param name="target_info">Information object, which should receive the time range and time values</param>
        void copy_time_information(vtkInformation* source_info, vtkInformation* target_info);

        /// <summary>
        /// Copy all time information (time range, time values, and current time step)
        /// </summary>
        /// <param name="source_alg">VTK algorithm object, which holds the current time step</param>
        /// <param name="target_alg">VTK algorithm object, which should receive the current time step</param>
        template <typename... alg_t>
        void copy_time_information(vtkDataSet* source_alg, vtkDataSet* target_alg, alg_t... target_algs);

        /// <summary>
        /// Copy all time information (time range, time values, and current time step)
        /// </summary>
        /// <param name="source_alg">VTK algorithm object, which holds the current time step</param>
        /// <param name="target_alg">VTK algorithm object, which should receive the current time step</param>
        void copy_time_information(vtkDataSet* source_alg, vtkDataSet* target_alg);
    }
}

#include "tpf_time.inl"
