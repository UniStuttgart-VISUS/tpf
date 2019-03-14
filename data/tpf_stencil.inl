#include "tpf_stencil.h"

#include "tpf_grid.h"

#include "../log/tpf_log.h"

#include "Eigen/Dense"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tpf
{
    namespace data
    {
        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline stencil<value_t, point_t, dimensions, rows, columns>::stencil(grid_type& grid, const extent_type& extent, behavior_t behavior, element_type constant)
            : data_grid(grid), extent(extent), behavior(behavior), constant(constant), temp_constant(constant) {}

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline stencil<value_t, point_t, dimensions, rows, columns>::stencil(grid_type& grid, const coords_type& location, std::size_t size, behavior_t behavior, element_type constant)
            : data_grid(grid), behavior(behavior), constant(constant), temp_constant(constant)
        {
            if (size % 2 != 1)
            {
                throw std::runtime_error(__tpf_error_message("Stencil size must be an odd number."));
            }
            if (behavior == behavior_t::GRADIENT && size != 3)
            {
                throw std::runtime_error(__tpf_error_message("Stencil size must be 3 for gradient behavior."));
            }

            this->extent.resize(dimensions);
            const std::size_t offset = (size - 1) / 2;

            for (std::size_t d = 0; d < dimensions; ++d)
            {
                this->extent[d].first = location[d] - offset;
                this->extent[d].second = location[d] + offset;
            }
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline std::size_t stencil<value_t, point_t, dimensions, rows, columns>::get_num_components() const
        {
            return rows * columns;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline std::size_t stencil<value_t, point_t, dimensions, rows, columns>::get_num_rows() const
        {
            return rows;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline std::size_t stencil<value_t, point_t, dimensions, rows, columns>::get_num_columns() const
        {
            return columns;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline std::size_t stencil<value_t, point_t, dimensions, rows, columns>::get_dimension() const
        {
            return dimensions;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline std::size_t stencil<value_t, point_t, dimensions, rows, columns>::get_num_elements() const
        {
            std::size_t num_elements = static_cast<std::size_t>(1);

            for (const auto& loc : this->location)
            {
                num_elements *= loc.second - loc.first + 1;
            }

            return num_elements;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline const typename stencil<value_t, point_t, dimensions, rows, columns>::extent_type& stencil<value_t, point_t, dimensions, rows, columns>::get_extent() const
        {
            return this->extent;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline void stencil<value_t, point_t, dimensions, rows, columns>::set_extent(const extent_type& extent)
        {
            this->extent = extent;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline void stencil<value_t, point_t, dimensions, rows, columns>::move(const Eigen::Matrix<long long, dimensions, 1> offset)
        {
            for (std::size_t d = 0; d < dimensions; ++d)
            {
                this->extent[d].first += offset[d];
                this->extent[d].second += offset[d];
            }
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline bool stencil<value_t, point_t, dimensions, rows, columns>::is_on_stencil(const coords_type& coords) const
        {
            bool on_stencil = true;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                on_stencil &= coords[i] >= 0LL && coords[i] <= this->extent[i].second - this->extent[i].first;
            }

            return on_stencil;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline const typename stencil<value_t, point_t, dimensions, rows, columns>::return_type stencil<value_t, point_t, dimensions, rows, columns>::operator()(const coords_type& coords) const
        {
            return const_cast<stencil<value_t, point_t, dimensions, rows, columns>*>(this)->access(coords);
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline typename stencil<value_t, point_t, dimensions, rows, columns>::return_type stencil<value_t, point_t, dimensions, rows, columns>::operator()(const coords_type& coords)
        {
            return access(coords);
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline Eigen::Matrix<point_t, dimensions, 1> stencil<value_t, point_t, dimensions, rows, columns>::get_cell_sizes(const coords_type& coords) const
        {
#ifdef __tpf_range_checks
            if (!is_on_stencil(coords))
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal stencil coordinates."));
            }
#endif

#ifdef __tpf_sanity_checks
            if (!this->data_grid.has_grid_information())
            {
                throw std::runtime_error(__tpf_error_message("Grid does not have cell sizes stored."));
            }
#endif

            // Calculate grid coordinates
            coords_type grid_coords = coords;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                grid_coords[i] += this->extent[i].first;
            }

            // Check: inside or outside the grid?
            if (this->data_grid.is_on_grid(grid_coords))
            {
                // Return element at the corresponding position inside the grid
                return this->data_grid.get_cell_sizes(grid_coords);
            }
            else
            {
                // Create temporary object with values according to the boundary behavior
                switch (this->behavior)
                {
                case MIRROR:
                    return this->data_grid.get_cell_sizes(mirror(grid_coords));

                    break;
                case CONSTANT:
                case REPEAT:
                case GRADIENT:
                default:
                    return this->data_grid.get_cell_sizes(clamp(grid_coords));
                }
            }
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline typename stencil<value_t, point_t, dimensions, rows, columns>::return_type stencil<value_t, point_t, dimensions, rows, columns>::access(const coords_type& coords)
        {
#ifdef __tpf_range_checks
            if (!is_on_stencil(coords))
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal stencil coordinates."));
            }
#endif

            // Calculate grid coordinates
            coords_type grid_coords = coords;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                grid_coords[i] += this->extent[i].first;
            }

            // Check: inside or outside the grid?
            if (this->data_grid.is_on_grid(grid_coords))
            {
                // Return element at the corresponding position inside the grid
                return this->data_grid(grid_coords);
            }
            else
            {
                // Create temporary object with values according to the boundary behavior
                switch (this->behavior)
                {
                case MIRROR:
                    return this->data_grid(mirror(grid_coords));

                    break;
                case REPEAT:
                    return this->data_grid(clamp(grid_coords));

                    break;
                case GRADIENT:
                    this->temp_constant = gradient(grid_coords);
                    return this->temp_constant;

                    break;
                case CONSTANT:
                default:
                    this->temp_constant = this->constant;
                    return this->temp_constant;
                }
            }
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline typename stencil<value_t, point_t, dimensions, rows, columns>::coords_type stencil<value_t, point_t, dimensions, rows, columns>::clamp(coords_type coords) const
        {
            for (std::size_t i = 0; i < dimensions; ++i)
            {
                coords[i] = std::max(coords[i], static_cast<long long>(this->data_grid.get_extent()[i].first));
                coords[i] = std::min(coords[i], static_cast<long long>(this->data_grid.get_extent()[i].second));
            }

            return coords;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline typename stencil<value_t, point_t, dimensions, rows, columns>::coords_type stencil<value_t, point_t, dimensions, rows, columns>::mirror(coords_type coords) const
        {
            for (std::size_t i = 0; i < dimensions; ++i)
            {
                if (coords[i] < static_cast<long long>(this->data_grid.get_extent()[i].first))
                {
                    coords[i] = 2 * static_cast<long long>(this->data_grid.get_extent()[i].first) - coords[i] - 1;
                }
                else if (coords[i] > static_cast<long long>(this->data_grid.get_extent()[i].second))
                {
                    coords[i] = 2 * static_cast<long long>(this->data_grid.get_extent()[i].second) - coords[i] + 1;
                }
            }

            return coords;
        }

        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline typename stencil<value_t, point_t, dimensions, rows, columns>::element_type stencil<value_t, point_t, dimensions, rows, columns>::gradient(coords_type coords) const
        {
            // This works only with stencil size 3, is checked in constructor.
            // Then if this function gets called, coords is within the first cell row outside the grid.
            const coords_type first = mirror(coords); // Mirror to get coords wihtin the first cell row in the grid.
            const coords_type second = 2 * first - coords; // Get the second cell row within the grid.

            point_t dist_first_second = static_cast<point_t>(1.0);
            point_t dist_first_coords = static_cast<point_t>(1.0);

            if (this->data_grid.has_grid_information())
            {
                dist_first_second = (this->data_grid.get_cell_coordinates(first) - this->data_grid.get_cell_coordinates(second)).norm();
                // data_grid.get_cell_coordinates() is only defined for coordinates within the grid.
                // According to get_cell_sizes() the cell size of the "coords" and "first" cell is the same.
                // Therefore the distance between them will match the cell size, but we need to distinguish
                // between side by side neighbors and diagonal neighbors.

                const auto first_size = this->data_grid.get_cell_sizes(first);
                dist_first_coords = static_cast<point_t>(0.0);

                for (std::size_t i = 0; i < dimensions; ++i)
                {
                    // Add cell length for all dimensions where first and coords differ. Square in loop for norm.
                    if (first[i] != coords[i])
                    {
                        dist_first_coords += first_size[i] * first_size[i];
                    }
                }
                dist_first_coords = std::sqrt(dist_first_coords);
            }

            const auto gradient = (this->data_grid(first) - this->data_grid(second)) / dist_first_second;
            const auto first_value = this->data_grid(first);

            return first_value + gradient * dist_first_coords;
        }
    }
}
