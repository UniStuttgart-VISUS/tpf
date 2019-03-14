#pragma once

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace math
    {
        /// <summary>
        /// Linear interpolation with given weight
        /// </summary>
        /// <template name="float_t">Floating point type for the weight</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="alpha">Interpolation weight in [0,1]</param>
        /// <param name="first">First value (0)</param>
        /// <param name="second">Second value (1)</param>
        /// <returns>Linear interpolated value</returns>
        template <typename float_t, typename value_t>
        value_t interpolate_linear(float_t alpha, const value_t& first, const value_t& second);

        /// <summary>
        /// Linear interpolation between given points
        /// </summary>
        /// <template name="float_t">Floating point type for points</template>
        /// <template name="dimension">Dimension of the points (allows linear interpolation in higher dimensions)</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="position">Position at which to interpolate</param>
        /// <param name="left">"Left" position (0)</param>
        /// <param name="right">"Right" position (1)</param>
        /// <param name="first">Value at the left position</param>
        /// <param name="second">Value at the right position</param>
        /// <returns>Linear interpolated value</returns>
        template <typename float_t, int dimension, typename value_t>
        value_t interpolate_linear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& left,
            const Eigen::Matrix<float_t, dimension, 1>& right, const value_t& first, const value_t& second);

        /// <summary>
        /// Bilinear interpolation with given weights
        /// </summary>
        /// <template name="float_t">Floating point type for the weight</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="alpha">Interpolation weight in [0,1] between first and second, and third and fourth</param>
        /// <param name="beta">Interpolation weight in [0,1] between linear interpolated values using alpha</param>
        /// <param name="first">First value (0,0)</param>
        /// <param name="second">Second value (1,0)</param>
        /// <param name="third">Third value (0,1)</param>
        /// <param name="fourth">Fourth value (1,1)</param>
        /// <returns>Bilinear interpolated value</returns>
        template <typename float_t, typename value_t>
        value_t interpolate_bilinear(float_t alpha, float_t beta,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth);

        /// <summary>
        /// Bilinear interpolation between given points
        /// </summary>
        /// <template name="float_t">Floating point type for points</template>
        /// <template name="dimension">Dimension of the points (allows linear interpolation in higher dimensions)</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="position">Position at which to interpolate</param>
        /// <param name="bottom_left">"Bottom left" position (0,0)</param>
        /// <param name="bottom_right">"Bottom right" position (1,0)</param>
        /// <param name="top_left">"Top left" position (0,1)</param>
        /// <param name="top_right">"Top right" position (1,1)</param>
        /// <param name="first">Value at the bottom left position</param>
        /// <param name="second">Value at the bottom right position</param>
        /// <param name="third">Value at the top left position</param>
        /// <param name="fourth">Value at the top right position</param>
        /// <returns>Bilinear interpolated value</returns>
        template <typename float_t, int dimension, typename value_t>
        value_t interpolate_bilinear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& bottom_right, const Eigen::Matrix<float_t, dimension, 1>& top_left, const Eigen::Matrix<float_t, dimension, 1>& top_right,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth);

        /// <summary>
        /// Trilinear interpolation with given weights
        /// </summary>
        /// <template name="float_t">Floating point type for the weight</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="alpha">Interpolation weight in [0,1] between first and second, and third and fourth</param>
        /// <param name="beta">Interpolation weight in [0,1] between linear interpolated values using alpha</param>
        /// <param name="gamma">Interpolation weight in [0,1] between linear interpolated values using beta</param>
        /// <param name="first">First value (0,0,0)</param>
        /// <param name="second">Second value (1,0,0)</param>
        /// <param name="third">Third value (0,1,0)</param>
        /// <param name="fourth">Fourth value (1,1,0)</param>
        /// <param name="fifth">Fifth value (0,0,1)</param>
        /// <param name="sixth">Sixth value (1,0,1)</param>
        /// <param name="seventh">Seventh value (0,1,1)</param>
        /// <param name="eighth">Eighth value (1,1,1)</param>
        /// <returns>Trilinear interpolated value</returns>
        template <typename float_t, typename value_t>
        value_t interpolate_trilinear(float_t alpha, float_t beta, float_t gamma,
            const value_t& first, const value_t& second, const value_t& third, const value_t& fourth,
            const value_t& fifth, const value_t& sixth, const value_t& seventh, const value_t& eighth);

        /// <summary>
        /// Trilinear interpolation between given points
        /// </summary>
        /// <template name="float_t">Floating point type for points</template>
        /// <template name="dimension">Dimension of the points (allows linear interpolation in higher dimensions)</template>
        /// <template name="value_t">Value type for the interpolation</template>
        /// <param name="position">Position at which to interpolate</param>
        /// <param name="front_bottom_left">"Front bottom left" position (0,0,0)</param>
        /// <param name="front_bottom_right">"Front bottom right" position (1,0,0)</param>
        /// <param name="front_top_left">"Front top left" position (0,1,0)</param>
        /// <param name="front_top_right">"Front top right" position (1,1,0)</param>
        /// <param name="back_bottom_left">"Back bottom left" position (0,0,1)</param>
        /// <param name="back_bottom_right">"Back bottom right" position (1,0,1)</param>
        /// <param name="back_top_left">"Back top left" position (0,1,1)</param>
        /// <param name="back_top_right">"Back top right" position (1,1,1)</param>
        /// <param name="first">Value at the front bottom left position</param>
        /// <param name="second">Value at the front bottom right position</param>
        /// <param name="third">Value at the front top left position</param>
        /// <param name="fourth">Value at the front top right position</param>
        /// <param name="fifth">Value at the back bottom left position</param>
        /// <param name="sixth">Value at the back bottom right position</param>
        /// <param name="seventh">Value at the back top left position</param>
        /// <param name="eighth">Value at the back top right position</param>
        /// <returns>Trilinear interpolated value</returns>
        template <typename float_t, int dimension, typename value_t>
        value_t interpolate_trilinear(const Eigen::Matrix<float_t, dimension, 1>& position, const Eigen::Matrix<float_t, dimension, 1>& front_bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& front_bottom_right, const Eigen::Matrix<float_t, dimension, 1>& front_top_left,
            const Eigen::Matrix<float_t, dimension, 1>& front_top_right, const Eigen::Matrix<float_t, dimension, 1>& back_bottom_left,
            const Eigen::Matrix<float_t, dimension, 1>& back_bottom_right, const Eigen::Matrix<float_t, dimension, 1>& back_top_left,
            const Eigen::Matrix<float_t, dimension, 1>& back_top_right, const value_t& first, const value_t& second, const value_t& third,
            const value_t& fourth, const value_t& fifth, const value_t& sixth, const value_t& seventh, const value_t& eighth);
    }
}

#include "tpf_interpolation.inl"
