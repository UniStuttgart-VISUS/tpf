#pragma once

#include "../traits/tpf_function.h"

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
        /// Policy for data structures supporting interpolation
        /// </summary>
        /// <template name="value_t">Value type, or return type of interpolation</template>
        /// <template name="point_t">Point type for interpolation position</template>
        template <typename value_t, typename point_t>
        class interpolatable
        {
        public:
            using interpolation_value_type = value_t;
            using interpolation_point_type = point_t;

            static constexpr bool value_type_interpolatable = has_operator_mul<float, value_t>::value && has_operator_add<value_t>::value;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate(const point_t& point) const = 0;
        };

        /// <summary>
        /// Policy for data structures supporting interpolation using Eigen::Matrix
        /// </summary>
        /// <template name="value_t">Value type of the Eigen::Matrix</template>
        /// <template name="rows">Number of rows of the Eigen::Matrix</template>
        /// <template name="columns">Number of columns of the Eigen::Matrix</template>
        /// <template name="point_t">Point type for interpolation position</template>
        template <typename value_t, int rows, int columns, typename point_t>
        class interpolatable<Eigen::Matrix<value_t, rows, columns>, point_t>
        {
        public:
            using interpolation_value_type = Eigen::Matrix<value_t, rows, columns>;
            using interpolation_point_type = point_t;

            static constexpr bool value_type_interpolatable = true;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual Eigen::Matrix<value_t, rows, columns> interpolate(const point_t& point) const = 0;
        };
    }
}