#pragma once

#include "Eigen/Dense"

namespace tpf
{
    namespace math
    {
        /// Convenience typedefs for matrices
        template <typename floatp_t> using mat2_t = Eigen::Matrix<floatp_t, 2, 2>;
        template <typename floatp_t> using mat3_t = Eigen::Matrix<floatp_t, 3, 3>;
        template <typename floatp_t> using mat4_t = Eigen::Matrix<floatp_t, 4, 4>;
    }
}