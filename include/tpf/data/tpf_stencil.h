#pragma once

#include "tpf_grid.h"

#include "Eigen/Dense"

#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Stencil from a rectilinear grid
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimensions">Number of spatial dimensions</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        template <typename value_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns = 1>
        class stencil
        {
        public:
            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(rows > 0, "Number of row components must be larger than zero");
            static_assert(columns > 0, "Number of column components must be larger than zero");

            /// Boundary behavior
            enum behavior_t
            {
                CONSTANT, MIRROR, REPEAT, GRADIENT
            };

            /// Used value and point type
            using value_type = typename std::remove_const<value_t>::type;
            using point_type = point_t;

            using grid_type = typename std::conditional<std::is_const<value_t>::value,
                const grid<value_type, point_t, dimensions, rows, columns>, grid<value_type, point_t, dimensions, rows, columns>>::type;

            using element_type = typename grid_type::element_type;
            using return_type = typename grid_type::return_type;

            /// Coordinate and extent type
            using coords_type = Eigen::Matrix<long long, dimensions, 1>;
            using extent_type = std::vector<std::pair<long long, long long>>;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="grid">Input grid</param>
            /// <param name="extent">Extent of the stencil</param>
            /// <param name="behavior">Border behavior</param>
            /// <param name="constant">Constant for constant border behavior</param>
            stencil(grid_type& grid, const extent_type& extent, behavior_t behavior = CONSTANT, element_type constant = static_cast<value_type>(0) * element_type());

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="grid">Input grid</param>
            /// <param name="location">Location of the stencil center</param>
            /// <param name="size">Stencil size per dimension</param>
            /// <param name="behavior">Border behavior</param>
            /// <param name="constant">Constant for constant border behavior</param>
            stencil(grid_type& grid, const coords_type& location, std::size_t size, behavior_t behavior = CONSTANT, element_type constant = static_cast<value_type>(0) * element_type());

            /// <summary>
            /// Get number of components
            /// </summary>
            /// <return>Number of components</return>
            std::size_t get_num_components() const;

            /// <summary>
            /// Get number of rows
            /// </summary>
            /// <return>Number of rows</return>
            std::size_t get_num_rows() const;

            /// <summary>
            /// Get number of columns
            /// </summary>
            /// <return>Number of columns</return>
            std::size_t get_num_columns() const;

            /// <summary>
            /// Get number of dimensions
            /// </summary>
            /// <return>Number of dimensions</return>
            std::size_t get_dimension() const;

            /// <summary>
            /// Get number of elements
            /// </summary>
            /// <return>Number of elements</return>
            std::size_t get_num_elements() const;

            /// <summary>
            /// Get extent of the stencil
            /// </summary>
            /// <returns>Extent of the stencil</returns>
            const extent_type& get_extent() const;

            /// <summary>
            /// Set extent of the stencil
            /// </summary>
            /// <param name="extent">Extent of the stencil</param>
            void set_extent(const extent_type& extent);

            /// <summary>
            /// Move stencil
            /// </summary>
            /// <param name="offset">Vector indicating the movement</param>
            void move(const Eigen::Matrix<long long, dimensions, 1> offset);

            /// <summary>
            /// Are the given coordinates on the stencil?
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>True if coordinates are on the stencil, false otherwise</return>
            bool is_on_stencil(const coords_type& coords) const;

            /// <summary>
            /// Are the given coordinates on the underlying grid?
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>True if (stencil) coordinates are on the underlying grid, false otherwise</return>
            bool is_on_grid(const coords_type& coords) const;

            /// <summary>
            /// Access data at the coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Data</return>
            const return_type operator()(const coords_type& coords) const;

            /// <summary>
            /// Access data at the coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Data</return>
            return_type operator()(const coords_type& coords);

            /// <summary>
            /// Get cell sizes at given coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Cell sizes at given coordinates</return>
            Eigen::Matrix<point_t, dimensions, 1> get_cell_sizes(const coords_type& coords) const;

        private:
            /// <summary>
            /// Access data at the coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <returns>Data</returns>
            return_type access(const coords_type& coords);

            /// <summary>
            /// Clamp the coordinates to the extent of the grid
            /// </summary>
            /// <param name="coords">Input coordinates</param>
            /// <returns>Clamped coordinates</returns>
            coords_type clamp(coords_type coords) const;

            /// <summary>
            /// Mirror the coordinates along the extent of the grid
            /// </summary>
            /// <param name="coords">Input coordinates</param>
            /// <returns>Mirrored coordinates</returns>
            coords_type mirror(coords_type coords) const;

            /// <summary>
            /// Calculate the value according to the gradient
            /// </summary>
            /// <param name="coords">Input coordinates</param>
            /// <returns>Calculate value</returns>
            element_type gradient(coords_type coords) const;

            /// Location and extent for each dimension
            extent_type extent;

            /// Grid
            grid_type& data_grid;

            /// Temporary value (for access to values outside the grid)
            mutable element_type temp_constant;

            /// Boundary behavior
            behavior_t behavior;
            element_type constant;
        };
    }
}

#include "tpf_stencil.inl"
