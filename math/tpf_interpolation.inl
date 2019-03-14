#include "tpf_interpolation.h"

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace math
    {
        template <typename float_t, typename value_t>
        inline value_t interpolate_linear(const float_t alpha, const value_t& first, const value_t& second)
        {
            return alpha * first + (static_cast<float_t>(1.0) - alpha) * second;
        }

        template <typename float_t, int dimension, typename value_t>
        inline value_t interpolate_linear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& left,
            const Eigen::Matrix<float_t, dimension, 1>& right, const value_t& first, const value_t& second)
        {
            static_assert(dimension > 0, "Dimension for linear interpolation must be at least 1");

            const auto vector_left_right = right - left;

            if (!(right - left).isZero())
            {
                const auto vector_left_position = position - left;

                const auto projected_position = left + (vector_left_position.dot(vector_left_right) / vector_left_right.dot(vector_left_right)) * vector_left_right;

                return interpolate_linear((projected_position - right).norm() / (left - right).norm(), first, second);
            }

            return first;
        }

        template <typename float_t, typename value_t>
        inline value_t interpolate_bilinear(const float_t alpha, const float_t beta,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth)
        {
            return interpolate_linear(beta, interpolate_linear(alpha, first, second), interpolate_linear(alpha, third, fourth));
        }

        template <typename float_t, int dimension, typename value_t>
        inline value_t interpolate_bilinear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& bottom_right, const Eigen::Matrix<float_t, dimension, 1>& top_left, const Eigen::Matrix<float_t, dimension, 1>& top_right,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth)
        {
            static_assert(dimension > 1, "Dimension for bilinear interpolation must be at least 2");

            return interpolate_linear(position, top_left, bottom_left,
                interpolate_linear(position, top_left, top_right, third, fourth),
                interpolate_linear(position, bottom_left, bottom_right, first, second));
        }

        template <typename float_t, typename value_t>
        inline value_t interpolate_trilinear(const float_t alpha, const float_t beta, const float_t gamma,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth,
            const value_t& fifth, const value_t& sixth, const value_t& seventh, const value_t& eighth)
        {
            return interpolate_linear(gamma, interpolate_bilinear(alpha, beta, first, second, third, fourth), interpolate_bilinear(alpha, beta, fifth, sixth, seventh, eighth));
        }

        template <typename float_t, int dimension, typename value_t>
        inline value_t interpolate_trilinear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& front_bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& front_bottom_right, const Eigen::Matrix<float_t, dimension, 1>& front_top_left,
            const Eigen::Matrix<float_t, dimension, 1>& front_top_right, const Eigen::Matrix<float_t, dimension, 1>& back_bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& back_bottom_right, const Eigen::Matrix<float_t, dimension, 1>& back_top_left,
            const Eigen::Matrix<float_t, dimension, 1>& back_top_right, const value_t& first, const value_t& second, const value_t& third,
            const value_t& fourth, const value_t& fifth, const value_t& sixth, const value_t& seventh, const value_t& eighth)
        {
            static_assert(dimension > 2, "Dimension for trilinear interpolation must be at least 3");

            return interpolate_linear(position, front_bottom_left, back_bottom_left,
                interpolate_bilinear(position, front_bottom_left, front_bottom_right, front_top_left, front_top_right, first, second, third, fourth),
                interpolate_bilinear(position, back_bottom_left, back_bottom_right, back_top_left, back_top_right, fifth, sixth, seventh, eighth));
        }
    }
}
