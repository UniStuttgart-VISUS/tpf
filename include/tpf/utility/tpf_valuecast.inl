#include "tpf_valuecast.h"

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace traits
    {
        template <typename new_t, typename original_t>
        constexpr inline new_t value_cast(original_t original)
        {
            return static_cast<new_t>(original);
        }

        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr inline Eigen::Matrix<new_value_t, rows, columns> value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
            typename std::enable_if<!std::is_same<original_value_t, new_value_t>::value>::type*)
        {
            return original.template cast<new_value_t>();
        }

        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr inline const Eigen::Matrix<new_value_t, rows, columns>& value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
            typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type*)
        {
            return original;
        }

        template <typename new_value_t, typename original_value_t, int rows, int columns>
        constexpr inline Eigen::Matrix<new_value_t, rows, columns> value_cast(Eigen::Matrix<original_value_t, rows, columns> original,
            typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type*)
        {
            return original;
        }
    }
}