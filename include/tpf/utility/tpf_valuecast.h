#pragma once

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace traits
    {
        /// <summary>
        /// Performs a cast using static_cast
        /// </summary>
        /// <template name="new_t">New type</template>
        /// <template name="original_t">Original type</template>
        /// <param name="original">Original value</param>
        /// <returns>Converted value</returns>
        template <typename new_t, typename original_t>
        constexpr new_t value_cast(original_t original);

        /// <summary>
        /// Performs a element-wise cast on vectors using static_cast
        /// </summary>
        /// <template name="new_value_t">New value type</template>
        /// <template name="original_value_t">Original value type</template>
        /// <param name="original">Original vector</param>
        /// <returns>Converted vector</returns>
        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr Eigen::Matrix<new_value_t, rows, columns> value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
            typename std::enable_if<!std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);

        /// <summary>
        /// Return the input unmodified
        /// </summary>
        /// <template name="new_value_t">New value type (same as original type)</template>
        /// <template name="original_value_t">Original value type</template>
        /// <param name="original">Original vector</param>
        /// <returns>Input vector</returns>
        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr const Eigen::Matrix<new_value_t, rows, columns>& value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
            typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);

        /// <summary>
        /// Return the input unmodified
        /// </summary>
        /// <template name="new_value_t">New value type (same as original type)</template>
        /// <template name="original_value_t">Original value type</template>
        /// <param name="original">Original vector</param>
        /// <returns>Input vector</returns>
        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr Eigen::Matrix<new_value_t, rows, columns> value_cast(Eigen::Matrix<original_value_t, rows, columns> original,
            typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);
    }
}

#include "tpf_valuecast.inl"