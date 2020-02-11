#pragma once

#include "vtkInformationVector.h"

namespace tpf
{
    namespace vtk
    {
        /// <summary>
        /// Set ghost level to required, or if higher to ghost level requested by filters downstream
        /// </summary>
        /// <param name="input_vector">Input information</param>
        /// <param name="output_vector">Output information</param>
        /// <param name="num_input_vectors">Number of input vectors</param>
        /// <param name="num_required_ghost_levels">Number of ghost levels required</param>
        /// <return>New number of ghost levels</return>
        std::size_t set_ghost_levels(vtkInformationVector **input_vector, vtkInformationVector *output_vector,
            const std::size_t num_input_vectors, const std::size_t num_required_ghost_levels);
    }
}

#include "tpf_ghost.inl"