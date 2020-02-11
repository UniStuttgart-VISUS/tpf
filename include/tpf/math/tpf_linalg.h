#pragma once

#include "Eigen/Dense"

#include <array>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace math
    {
        /// <summary>
        /// Struct for storing eigenpairs
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="dimension">Vector dimension</template>
        template <typename floatp_t, std::size_t dimension>
        struct eigenpairs
        {
            static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

            std::array<Eigen::Matrix<floatp_t, dimension, 1>, dimension> eigenvectors;
            std::array<floatp_t, dimension> eigenvalues;
        };

        /// <summary>
        /// Calculate gradient using a weighted regression based method
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimension">Vector dimension</template>
        /// <param name="vertices">Input vertices</param>
        /// <param name="values">Values corresponding to the input vertices</param>
        /// <param name="inverse_distance_weighting">Use inverse distance weighting?</param>
        /// <return>Jacobian matrix</return>
        template <typename value_t, typename point_t, std::size_t dimension>
        Eigen::Matrix<value_t, dimension, 1> calculate_gradient(const std::vector<Eigen::Matrix<point_t, dimension, 1>>& vertices,
            const std::vector<value_t>& values, const bool inverse_distance_weighting = true);

        /// <summary>
        /// Calculate Jacobian matrix using a weighted regression based method
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimension">Vector dimension</template>
        /// <param name="vertices">Input vertices</param>
        /// <param name="values">Values corresponding to the input vertices</param>
        /// <param name="inverse_distance_weighting">Use inverse distance weighting?</param>
        /// <return>Jacobian matrix</return>
        template <typename value_t, typename point_t, std::size_t dimension>
        Eigen::Matrix<value_t, dimension, dimension> calculate_jacobian(const std::vector<Eigen::Matrix<point_t, dimension, 1>>& vertices,
            const std::vector<Eigen::Matrix<value_t, dimension, 1>>& values, const bool inverse_distance_weighting = true);

        /// <summary>
        /// Extract symmetric part of a matrix
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="rows_and_columns">Number of matrix rows and columns</template>
        /// <param name="matrix">Input matrix</param>
        /// <return>Matrix containing only the symmetric part</return>
        template <typename floatp_t, std::size_t rows_and_columns>
        Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> extract_symmetric_part(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix);

        /// <summary>
        /// Extract asymmetric part of a matrix
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="rows_and_columns">Number of matrix rows and columns</template>
        /// <param name="matrix">Input matrix</param>
        /// <return>Matrix containing only the asymmetric part</return>
        template <typename floatp_t, std::size_t rows_and_columns>
        Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> extract_asymmetric_part(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix);

        /// <summary>
        /// Calculate eigenvalues and eigenvectors of a matrix
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="rows_and_columns">Number of matrix rows and columns</template>
        /// <param name="matrix">Input matrix</param>
        /// <return>Calculated eigenvalues and eigenvectors</return>
        template <typename floatp_t, std::size_t rows_and_columns>
        eigenpairs<floatp_t, rows_and_columns> calculate_eigenpair(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix);
    }
}

#include "tpf_linalg.inl"