#pragma once

#include "tpf_array.h"
#include "tpf_data_information.h"
#include "tpf_grid_information.h"

#include "../policies/tpf_interpolatable.h"

#include "../utility/tpf_optional.h"

#include "Eigen/Dense"

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Rectilinear grid
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimensions">Number of spatial dimensions</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        template <typename value_t, typename point_t, int dimensions = 3, int rows = 1, int columns = 1>
        class grid :
            public policies::interpolatable<typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type, Eigen::Matrix<point_t, dimensions, 1>>
        {
        public:
            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(rows > 0, "Number of row components must be larger than zero");
            static_assert(columns > 0, "Number of column components must be larger than zero");

            using interpolation_policy = policies::interpolatable<typename std::conditional<(rows * columns) != 1,
                Eigen::Matrix<value_t, rows, columns>, value_t>::type, Eigen::Matrix<point_t, dimensions, 1>>;

            /// Used value and point type
            using value_type = value_t;
            using point_type = point_t;

            using element_type = typename array<value_t, rows, columns>::element_type;
            using return_type = typename array<value_t, rows, columns>::return_type;

            /// Coordinate type
            using coords_t = Eigen::Matrix<long long, dimensions, 1>;

            /// <summary>
            /// Create a grid of empty size
            /// </summary>
            /// <remarks>Careful! Only get_size() has defined behavior</remarks>
            grid();

            /// <summary>
            /// Create grid with given extents and number of components
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="extent">Extent</param>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <param name="node_coordinates">Node coordinates</param>
            /// <param name="cell_sizes">Cell sizes</param>
            grid(std::string name, extent_t extent,
                typename grid_information<point_t>::array_type cell_coordinates = typename grid_information<point_t>::array_type(),
                typename grid_information<point_t>::array_type node_coordinates = typename grid_information<point_t>::array_type(),
                typename grid_information<point_t>::array_type cell_sizes = typename grid_information<point_t>::array_type());

            /// <summary>
            /// Create grid with given extents and number of components
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="extent">Extent</param>
            /// <param name="data_array">Data</param>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <param name="node_coordinates">Node coordinates</param>
            /// <param name="cell_sizes">Cell sizes</param>
            grid(std::string name, extent_t extent, std::vector<value_t> data_array,
                typename grid_information<point_t>::array_type cell_coordinates = typename grid_information<point_t>::array_type(),
                typename grid_information<point_t>::array_type node_coordinates = typename grid_information<point_t>::array_type(),
                typename grid_information<point_t>::array_type cell_sizes = typename grid_information<point_t>::array_type());

            /// <summary>
            /// Default copy/move constructor/assignment
            /// </summary>
            grid(const grid&) = default;
            grid(grid&&) noexcept = default;
            grid& operator=(const grid&) = default;
            grid& operator=(grid&&) noexcept = default;

            /// <summary>
            /// Initialize grid with given value
            /// </summary>
            /// <param name="value">Value for every entry</param>
            void initialize(const value_t& value) noexcept;

            /// <summary>
            /// Initialize grid with given vector
            /// </summary>
            /// <param name="value">Vector or matrix for every entry</param>
            template <int _enable_rows = rows, int _enable_columns = columns>
            void initialize(const element_type& value, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type* = nullptr) noexcept;

            /// <summary>
            /// Are the given coordinates on the grid?
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>True if coordinates are on the grid, false otherwise</return>
            bool is_on_grid(const coords_t& coords) const noexcept;

            /// <summary>
            /// Are the given coordinates inside (and not on the ghost levels)?
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="num_ghost_levels">Number of ghost levels</param>
            /// <return>True if coordinates are inside, false otherwise</return>
            bool is_local(const coords_t& coords, std::size_t num_ghost_levels) const noexcept;

            /// <summary>
            /// Are the given coordinates on a ghost level?
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="num_ghost_levels">Number of ghost levels</param>
            /// <return>True if coordinates are on a ghost level, false otherwise</return>
            bool is_ghost(const coords_t& coords, std::size_t num_ghost_levels) const noexcept;

            /// <summary>
            /// Find cell for a given vertex
            /// </summary>
            /// <param name="vertex">Vertex</param>
            /// <return>[valid, cell coordinates]</return>
            utility::optional<coords_t> find_cell(const Eigen::Matrix<point_t, dimensions, 1>& vertex) const;

            /// <summary>
            /// Check if the vertex is in the given cell
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="vertex">Vertex</param>
            /// <return>True if vertex is in the given cell, false otherwise</return>
            bool is_in_cell(const coords_t& coords, const Eigen::Matrix<point_t, dimensions, 1>& vertex) const;

            /// <summary>
            /// Find staggered cell for a given vertex
            /// </summary>
            /// <param name="vertex">Vertex</param>
            /// <return>[valid, cell coordinates]</return>
            utility::optional<coords_t> find_staggered_cell(const Eigen::Matrix<point_t, dimensions, 1>& vertex) const;

            /// <summary>
            /// Check if the vertex is in the given staggered cell
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="vertex">Vertex</param>
            /// <return>True if vertex is in the given cell, false otherwise</return>
            bool is_in_staggered_cell(const coords_t& coords, const Eigen::Matrix<point_t, dimensions, 1>& vertex) const;

            /// <summary>
            /// Find cell for a given polygon, evaluating the position at its barycenter
            /// </summary>
            /// <param name="vertices">Vertices of the polygon</param>
            /// <return>[valid, cell coordinates]</return>
            utility::optional<coords_t> find_cell(const std::vector<Eigen::Matrix<point_t, dimensions, 1>>& vertices) const;

            /// <summary>
            /// Get name
            /// </summary>
            /// <return>Name</return>
            const std::string& get_name() const noexcept;

            /// <summary>
            /// Set name
            /// </summary>
            /// <param name="name">Name</param>
            void set_name(std::string name) noexcept;

            /// <summary>
            /// Get extent
            /// </summary>
            /// <return>Extent</return>
            const extent_t& get_extent() const noexcept;

            /// <summary>
            /// Set grid information
            /// </summary>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <param name="node_coordinates">Node coordinates</param>
            /// <param name="cell_sizes">Cell sizes</param>
            void set_grid_information(typename grid_information<point_t>::array_type cell_coordinates,
                typename grid_information<point_t>::array_type node_coordinates, typename grid_information<point_t>::array_type cell_sizes);

            /// <summary>
            /// Get grid information
            /// </summary>
            /// <returns>Grid information</returns>
            grid_information<point_t> get_grid_information() const;

            /// <summary>
            /// Get cell coordinates
            /// </summary>
            /// <return>Cell coordinates</return>
            const typename grid_information<point_t>::array_type& get_cell_coordinates() const;

            /// <summary>
            /// Get cell coordinates at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Cell coordinates at given coordinates and for given dimension</return>
            const point_t& get_cell_coordinates(const coords_t& coords, std::size_t dimension) const;

            /// <summary>
            /// Get cell coordinates at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Cell coordinates at given coordinates and for given dimension</return>
            const point_t& get_cell_coordinates(std::size_t coords, std::size_t dimension) const;

            /// <summary>
            /// Get cell coordinates at given coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Cell coordinates at given coordinates</return>
            Eigen::Matrix<point_t, dimensions, 1> get_cell_coordinates(const coords_t& coords) const;

            /// <summary>
            /// Get node coordinates
            /// </summary>
            /// <return>Node coordinates</return>
            const typename grid_information<point_t>::array_type& get_node_coordinates() const;

            /// <summary>
            /// Get node coordinates at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Node coordinates at given coordinates and for given dimension</return>
            const point_t& get_node_coordinates(const coords_t& coords, std::size_t dimension) const;

            /// <summary>
            /// Get node coordinates at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Node coordinates at given coordinates and for given dimension</return>
            const point_t& get_node_coordinates(std::size_t coords, std::size_t dimension) const;

            /// <summary>
            /// Get node coordinates at given coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Node coordinates at given coordinates</return>
            Eigen::Matrix<point_t, dimensions, 1> get_node_coordinates(const coords_t& coords) const;

            /// <summary>
            /// Get cell sizes
            /// </summary>
            /// <return>Cell sizes</return>
            const typename grid_information<point_t>::array_type& get_cell_sizes() const;

            /// <summary>
            /// Get cell size at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Cell size at given coordinates and for given dimension</return>
            const point_t& get_cell_sizes(const coords_t& coords, std::size_t dimension) const;

            /// <summary>
            /// Get cell sizes at given coordinates and for given dimension
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="dimension">dimension</param>
            /// <return>Cell sizes at given coordinates and for given dimension</return>
            const point_t& get_cell_sizes(std::size_t coords, std::size_t dimension) const;

            /// <summary>
            /// Get cell sizes at given coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Cell sizes at given coordinates</return>
            Eigen::Matrix<point_t, dimensions, 1> get_cell_sizes(const coords_t& coords) const;

            /// <summary>
            /// Is there grid information?
            /// </summary>
            /// <returns>Grid information available?</returns>
            bool has_grid_information() const noexcept;

            /// <summary>
            /// Get number of components
            /// </summary>
            /// <return>Number of components</return>
            constexpr std::size_t get_num_components() const noexcept;

            /// <summary>
            /// Get number of rows
            /// </summary>
            /// <return>Number of rows</return>
            constexpr std::size_t get_num_rows() const noexcept;

            /// <summary>
            /// Get number of columns
            /// </summary>
            /// <return>Number of columns</return>
            constexpr std::size_t get_num_columns() const noexcept;

            /// <summary>
            /// Get dimension
            /// </summary>
            /// <return>dimension</return>
            constexpr std::size_t get_dimension() const noexcept;

            /// <summary>
            /// Get number of elements
            /// </summary>
            /// <return>Number of elements</return>
            std::size_t get_num_elements() const noexcept;

            /// <summary>
            /// Get array size
            /// </summary>
            /// <return>Array size</return>
            std::size_t get_size() const noexcept;

            /// <summary>
            /// Access data
            /// </summary>
            /// <return>Data</return>
            const std::vector<value_t>& get_data() const noexcept;

            /// <summary>
            /// Access data
            /// </summary>
            /// <return>Data</return>
            std::vector<value_t>& get_data() noexcept;

            /// <summary>
            /// Set data
            /// </summary>
            /// <param name="data_array">Data array</param>
            void set_data(std::vector<value_t> data_array);

            /// <summary>
            /// Access array
            /// </summary>
            /// <return>Array</return>
            const array<value_t, rows, columns>& get_array() const noexcept;

            /// <summary>
            /// Access array
            /// </summary>
            /// <return>Array</return>
            array<value_t, rows, columns>& get_array() noexcept;

            /// <summary>
            /// Access data at the coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Vector</return>
            const return_type operator()(const coords_t& coords) const;

            /// <summary>
            /// Access data at the coordinates
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Vector</return>
            return_type operator()(const coords_t& coords);

            /// <summary>
            /// Copy the structure of this grid to create a new one
            /// </summary>
            /// <template name="new_value_t">Alternative value type of the new grid</template>
            /// <template name="new_rows">Alternative number of row components of the new grid</template>
            /// <template name="new_columns">Alternative number of column components of the new grid</template>
            /// <param name="name">Name of the new grid</param>
            /// <returns>New grid</returns>
            template <typename new_value_t = value_t, int new_rows = rows, int new_columns = columns>
            grid<new_value_t, point_t, dimensions, new_rows, new_columns> copy_structure(std::string name) const;

            /// <summary>
            /// Extract a subgrid, given by the extent
            /// </summary>
            /// <template name="new_value_t">Alternative value type for the extracted subgrid</template>
            /// <template name="new_point_t">Alternative point type for the extracted subgrid</template>
            /// <param name="subextent">Extent of subgrid</param>
            /// <returns>Subgrid</returns>
            template <typename new_value_t = value_t, typename new_point_t = point_t>
            grid<new_value_t, new_point_t, dimensions, rows, columns> extract_subgrid(const extent_t& subextent) const;

            /// <summary>
            /// Extract subgrid information for the given extent
            /// </summary>
            /// <template name="new_value_t">Alternative point type for the extracted subgrid</template>
            /// <param name="subextent">Extent of subgrid</param>
            /// <returns>Subgrid information</returns>
            template <typename new_point_t = point_t>
            grid_information<new_point_t> extract_subgrid_information(const extent_t& subextent) const;

            /// <summary>
            /// Copy data from input grid into the extent specified
            /// </summary>
            /// <template name="new_value_t">Alternative value type for the input subgrid</template>
            /// <param name="input">Input grid</param>
            /// <param name="subextent">Extent of this grid to be replaced</param>
            /// <param name="fill_grid_info">Fill also grid information</param>
            template <typename new_value_t = value_t>
            void fill_subgrid(const grid<new_value_t, point_t, dimensions, rows, columns>& input, bool fill_grid_info = false);

            /// <summary>
            /// Copy data from input grid into the extent specified
            /// </summary>
            /// <template name="new_value_t">Alternative value type for the input subgrid</template>
            /// <param name="input">Input grid</param>
            /// <param name="subextent">Extent of this grid to be replaced</param>
            /// <param name="fill_grid_info">Fill also grid information</param>
            template <typename new_value_t = value_t>
            void fill_subgrid(const grid<new_value_t, point_t, dimensions, rows, columns>& input, const extent_t& subextent, bool fill_grid_info = false);

            /// <summary>
            /// Interpolate value at a position within the grid
            /// If the underlying value type does not support interpolation, the nearest neighbor is returned.
            /// </summary>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            virtual typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type
                interpolate(const Eigen::Matrix<point_t, dimensions, 1>& point) const override;
        private:
            /// <summary>
            /// Calculate the index for given coordinates and component
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <return>Index</return>
            std::size_t calculate_index(const coords_t& coords) const;

            /// <summary>
            /// Interpolate value at a position within the grid
            /// </summary>
            /// <template name="active">Indicates whether interpolation is supported</template>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            template <bool active = interpolation_policy::value_type_interpolatable>
            typename std::enable_if<active, typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type>::type
                interpolate_impl(const Eigen::Matrix<point_t, dimensions, 1>& point) const;

            /// <summary>
            /// Return nearest neighbor for a position within the grid
            /// </summary>
            /// <template name="active">Indicates whether interpolation is supported</template>
            /// <param name="point">Point at which to interpolate</param>
            /// <returns>Interpolated value</returns>
            template <bool active = interpolation_policy::value_type_interpolatable>
            typename std::enable_if<!active, typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type>::type
                interpolate_impl(const Eigen::Matrix<point_t, dimensions, 1>& point) const;

            /// Extent for each dimension
            extent_t extent;

            /// Cell coordinates
            typename grid_information<point_t>::array_type cell_coordinates;

            /// Node coordinates
            typename grid_information<point_t>::array_type node_coordinates;

            /// Cell sizes
            typename grid_information<point_t>::array_type cell_sizes;

            /// Is there grid information?
            bool grid_information_avail;

            /// Data
            array<value_t, rows, columns> data;
        };
    }
}

#include "tpf_grid.inl"