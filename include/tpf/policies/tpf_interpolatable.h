#pragma once

#include "../traits/tpf_function.h"

#include <vector>

namespace tpf
{
    namespace policies
    {
        namespace
        {
            __tpf_has_operator(operator*, *, mul);
            __tpf_has_operator(operator+, +, add);
        }

        /// <summary>
        /// Base class of the policy for data structures supporting interpolation
        /// </summary>
        /// <template name="point_t">Point type for interpolation position</template>
        template <typename point_t>
        class interpolatable_base
        {
        public:
            /// <summary>
            /// Return number of components of the interpolated values
            /// </summary>
            /// <returns>Number of components of interpolated values</returns>
            virtual std::size_t get_num_interpolated_components() const = 0;

            /// <summary>
            /// Interpolate at a certain point and return interpolated value as a vector of doubles
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual std::vector<double> interpolate_dynamic(const point_t& point) const = 0;
        };

        /// <summary>
        /// Policy for data structures supporting interpolation
        /// </summary>
        /// <template name="value_t">Value type, or return type of interpolation</template>
        /// <template name="point_t">Point type for interpolation position</template>
        template <typename value_t, typename point_t>
        class interpolatable : public interpolatable_base<point_t>
        {
        public:
            using interpolation_value_type = value_t;
            using interpolation_point_type = point_t;

            static constexpr bool value_type_interpolatable = has_operator_mul<float, value_t>::value && has_operator_add<value_t>::value;

            /// <summary>
            /// Return number of components of the interpolated values
            /// </summary>
            /// <returns>Number of components of interpolated values</returns>
            virtual std::size_t get_num_interpolated_components() const override;

            /// <summary>
            /// Interpolate at a certain point and return interpolated value
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate(const point_t& point) const = 0;

            /// <summary>
            /// Interpolate at a certain point and return interpolated value as a vector of doubles
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual std::vector<double> interpolate_dynamic(const point_t& point) const override;
        };

        /// <summary>
        /// Policy for data structures supporting interpolation using Eigen::Matrix
        /// </summary>
        /// <template name="value_t">Value type of the Eigen::Matrix</template>
        /// <template name="rows">Number of rows of the Eigen::Matrix</template>
        /// <template name="columns">Number of columns of the Eigen::Matrix</template>
        /// <template name="point_t">Point type for interpolation position</template>
        template <typename value_t, int rows, int columns, typename point_t>
        class interpolatable<Eigen::Matrix<value_t, rows, columns>, point_t> : public interpolatable_base<point_t>
        {
        public:
            using interpolation_value_type = Eigen::Matrix<value_t, rows, columns>;
            using interpolation_point_type = point_t;

            static constexpr bool value_type_interpolatable = true;

            /// <summary>
            /// Return number of components of the interpolated values
            /// </summary>
            /// <returns>Number of components of interpolated values</returns>
            virtual std::size_t get_num_interpolated_components() const override;

            /// <summary>
            /// Interpolate at a certain point and return interpolated value
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual Eigen::Matrix<value_t, rows, columns> interpolate(const point_t& point) const = 0;

            /// <summary>
            /// Interpolate at a certain point and return interpolated value as a vector of doubles
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual std::vector<double> interpolate_dynamic(const point_t& point) const override;
        };
    }
}

#include "tpf_interpolatable.inl"
