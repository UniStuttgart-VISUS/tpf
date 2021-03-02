#pragma once

#include "tpf_interpolatable.h"

#include <vector>

namespace tpf
{
    namespace policies
    {
        template <typename value_t, typename point_t>
        inline std::size_t interpolatable<value_t, point_t>::get_num_interpolated_components() const
        {
            return 1;
        }

        template <typename value_t, typename point_t>
        inline std::vector<double> interpolatable<value_t, point_t>::interpolate_dynamic(const point_t& point) const
        {
            return std::vector<double>{ static_cast<double>(interpolate(point)) };
        }

        template <typename value_t, int rows, int columns, typename point_t>
        inline std::size_t interpolatable<Eigen::Matrix<value_t, rows, columns>, point_t>::get_num_interpolated_components() const
        {
            return rows * columns;
        }

        template <typename value_t, int rows, int columns, typename point_t>
        inline std::vector<double> interpolatable<Eigen::Matrix<value_t, rows, columns>, point_t>::interpolate_dynamic(const point_t& point) const
        {
            const auto interpolated = interpolate(point);

            std::vector<double> value(rows * columns);

            for (std::size_t j = 0; j < columns; ++j)
            {
                for (std::size_t i = 0; i < rows; ++i)
                {
                    value[j * rows + i] = interpolated(i, j);
                }
            }

            return value;
        }
    }
}
