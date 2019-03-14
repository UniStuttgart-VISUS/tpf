#pragma once

#include "tpf_math_defines.h"

#include "Eigen/Dense"

#include <array>
#include <utility>

namespace tpf
{
    namespace math
    {
        /// Convenience typedefs for vectors
        template <typename floatp_t> using vec2_t = Eigen::Matrix<floatp_t, 2, 1>;
        template <typename floatp_t> using vec3_t = Eigen::Matrix<floatp_t, 3, 1>;
        template <typename floatp_t> using vec4_t = Eigen::Matrix<floatp_t, 4, 1>;

        /// <summary>
        /// Get two normalized vectors with directions orthogonal to the given vector
        /// </summary>
        /// <template name="floatp_t">Floating point type for vector values</template>
        /// <param name="vector">Vector to which orthogonal ones should be found</param>
        /// <return>Pair of orthogonal vectors</return>
        template <typename floatp_t>
        std::pair<Eigen::Matrix<floatp_t, 3, 1>, Eigen::Matrix<floatp_t, 3, 1>> orthonormal(const Eigen::Matrix<floatp_t, 3, 1>& vector, floatp_t epsilon = default_epsilon<floatp_t>::value);

        /// <summary>
        /// Calculate angle between two vectors
        /// </summary>
        /// <template name="floatp_t">Floating point type for vector values</template>
        /// <param name="first">First vector</param>
        /// <param name="second">Second vector</param>
        /// <return>Angle spanned by the two vectors</return>
        template <typename floatp_t>
        floatp_t calculate_angle(const Eigen::Matrix<floatp_t, 3, 1>& first, const Eigen::Matrix<floatp_t, 3, 1>& second, floatp_t epsilon = default_epsilon<floatp_t>::value);

        /// <summary>
        /// Find directions in order of similarity to the vector
        /// </summary>
        /// <template name="floatp_t">Floating point type for vector values</template>
        /// <param name="vector">Vector to which fitting directions should be found</param>
        /// <return>Three best-fitting directions</return>
        template <typename floatp_t, int components>
        std::array<Eigen::Matrix<floatp_t, components, 1>, components> get_fitting_directions(const Eigen::Matrix<floatp_t, components, 1>& vector);

        /// <summary>
        /// Find direction most similar to the vector
        /// </summary>
        /// <template name="floatp_t">Floating point type for vector values</template>
        /// <param name="vector">Vector to which fitting directions should be found</param>
        /// <return>Best-fitting direction</return>
        template <typename floatp_t, int components>
        Eigen::Matrix<floatp_t, components, 1> get_fitting_direction(const Eigen::Matrix<floatp_t, 3, 1>& vector);

        /// <summary>
        /// Test two vectors for equality within an error margin per component
        /// </summary>
        /// <param name="lhs">First vector for comparison</param>
        /// <param name="rhs">Second vector for comparison</param>
        /// <param name="epsilon">Error margin allowed</param>
        /// <returns>True if equal within an error margin; false otherwise</returns>
        template <typename floatp_t, int components>
        bool equals(const Eigen::Matrix<floatp_t, components, 1>& lhs, const Eigen::Matrix<floatp_t, components, 1>& rhs, floatp_t epsilon = default_epsilon<floatp_t>::value);
    }
}

#include "tpf_vector.inl"
