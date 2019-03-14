#include "tpf_grid.h"

#include "tpf_data_information.h"
#include "tpf_grid_information.h"

#include "../algorithm/tpf_loop.h"

#include "../log/tpf_log.h"

#include "../math/tpf_interpolation.h"

#include "../utility/tpf_optional.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline grid<value_t, point_t, dimensions, rows, columns>::grid() : data("EMPTY"), grid_information_avail(false)
        { }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline grid<value_t, point_t, dimensions, rows, columns>::grid(std::string name, extent_t extent, typename grid_information<point_t>::array_type cell_coordinates,
            typename grid_information<point_t>::array_type node_coordinates, typename grid_information<point_t>::array_type cell_sizes)
                : data(std::move(name)), grid_information_avail(false)
        {
            if (extent.size() != dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Extent dimension doesn't match the number of dimensions specified."));
            }

            std::swap(this->extent, extent);

            std::size_t num_items = 1;

            for (const auto& ext : this->extent)
            {
                num_items *= ext.second - ext.first + 1;
            }

            this->data.resize(num_items);

            set_grid_information(std::move(cell_coordinates), std::move(node_coordinates), std::move(cell_sizes));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline grid<value_t, point_t, dimensions, rows, columns>::grid(std::string name, extent_t extent,
            std::vector<value_t> data_array, typename grid_information<point_t>::array_type cell_coordinates, typename grid_information<point_t>::array_type node_coordinates,
            typename grid_information<point_t>::array_type cell_sizes) : data(std::move(name), std::move(data_array)), grid_information_avail(false)
        {
            if (extent.size() != dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Extent dimension doesn't match the number of dimensions specified."));
            }

            std::swap(this->extent, extent);

            std::size_t num_items = 1;

            for (const auto& ext : this->extent)
            {
                num_items *= ext.second - ext.first + 1;
            }

            if (num_items != this->data.get_num_elements())
            {
                throw std::runtime_error(__tpf_error_message("Data array size doesn't match the extent specified."));
            }

            set_grid_information(std::move(cell_coordinates), std::move(node_coordinates), std::move(cell_sizes));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline void grid<value_t, point_t, dimensions, rows, columns>::initialize(const value_t& value) noexcept
        {
            this->data.initialize(value);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline void grid<value_t, point_t, dimensions, rows, columns>::initialize(const element_type& value, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type*) noexcept
        {
            this->data.initialize(value);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline bool grid<value_t, point_t, dimensions, rows, columns>::is_on_grid(const coords_t& coords) const noexcept
        {
            bool on_grid = true;

            for (int i = 0; i < dimensions; ++i)
            {
                on_grid &= coords[i] >= 0LL && static_cast<std::size_t>(coords[i]) >= this->extent[i].first && static_cast<std::size_t>(coords[i]) <= this->extent[i].second;
            }

            return on_grid;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline bool grid<value_t, point_t, dimensions, rows, columns>::is_local(const coords_t& coords, const std::size_t num_ghost_levels) const noexcept
        {
            bool inside = true;

            for (int i = 0; i < dimensions; ++i)
            {
                inside &= coords[i] >= 0LL && static_cast<std::size_t>(coords[i]) >= this->extent[i].first + num_ghost_levels
                    && static_cast<std::size_t>(coords[i]) <= this->extent[i].second - num_ghost_levels;
            }

            return inside;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline bool grid<value_t, point_t, dimensions, rows, columns>::is_ghost(const coords_t& coords, const std::size_t num_ghost_levels) const noexcept
        {
            return is_on_grid(coords) && !is_local(coords, num_ghost_levels);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline utility::optional<typename grid<value_t, point_t, dimensions, rows, columns>::coords_t> grid<value_t, point_t, dimensions, rows, columns>::find_cell(
            const Eigen::Matrix<point_t, dimensions, 1>& vertex) const
        {
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Grid information has not been set."));
            }

            // Find cell
            coords_t coords;
            std::vector<bool> valid(dimensions);

            for (std::size_t c = 0; c < dimensions; ++c)
            {
                valid[c] = false;

                for (std::size_t i = 0; i < this->node_coordinates[c].size() - 1; ++i)
                {
                    if (this->node_coordinates[c][i] <= vertex[c] && this->node_coordinates[c][i + 1] >= vertex[c])
                    {
                        coords[c] = i + this->extent[c].first;
                        valid[c] = true;
                    }
                }
            }

            bool all_valid = true;

            for (std::size_t c = 0; c < dimensions; ++c)
            {
                all_valid &= valid[c];
            }

            if (all_valid)
            {
                return coords;
            }
            else
            {
                return utility::nullopt;
            }
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline utility::optional<typename grid<value_t, point_t, dimensions, rows, columns>::coords_t> grid<value_t, point_t, dimensions, rows, columns>::find_cell
            (const std::vector<Eigen::Matrix<point_t, dimensions, 1>>& vertices) const
        {
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Grid information has not been set."));
            }

            // Get point inside of the cell (if possible)
            Eigen::Matrix<point_t, dimensions, 1> center;

            if (vertices.size() == 1)
            {
                center = vertices[0];
            }
            else if (vertices.size() == 2)
            {
                center = vertices[0] + static_cast<point_t>(0.5L) * (vertices[1] - vertices[0]);
            }
            else if (vertices.size() == 3)
            {
                center = static_cast<point_t>(1.0L / 3.0L) * (vertices[0] + vertices[1] + vertices[2]);
            }
            else
            {
                center = vertices[0] + static_cast<point_t>(0.5L) * (vertices[vertices.size() / 2] - vertices[0]);
            }

            // Find cell
            return find_cell(center);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const std::string& grid<value_t, point_t, dimensions, rows, columns>::get_name() const noexcept
        {
            return this->data.get_name();
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline void grid<value_t, point_t, dimensions, rows, columns>::set_name(std::string name) noexcept
        {
            this->data.set_name(std::move(name));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const extent_t& grid<value_t, point_t, dimensions, rows, columns>::get_extent() const noexcept
        {
            return this->extent;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline void grid<value_t, point_t, dimensions, rows, columns>::set_grid_information(typename grid_information<point_t>::array_type cell_coordinates,
            typename grid_information<point_t>::array_type node_coordinates, typename grid_information<point_t>::array_type cell_sizes)
        {
            if (cell_coordinates.size() != 0 || node_coordinates.size() != 0 || cell_sizes.size() != 0)
            {
                if (cell_coordinates.size() != dimensions || node_coordinates.size() != dimensions || cell_sizes.size() != dimensions)
                {
                    throw std::runtime_error(__tpf_error_message("Coordinate dimension doesn't match the number of dimensions specified."));
                }

                for (std::size_t i = 0; i < dimensions; ++i)
                {
                    std::size_t num_elements = extent[i].second - extent[i].first + 1;

                    if (cell_coordinates[i].size() != num_elements || node_coordinates[i].size() != num_elements + 1 || cell_sizes[i].size() != num_elements)
                    {
                        throw std::runtime_error(__tpf_error_message("Coordinate extent doesn't match the extent specified."));
                    }
                }

                this->grid_information_avail = true;
            }

            std::swap(this->cell_coordinates, cell_coordinates);
            std::swap(this->node_coordinates, node_coordinates);
            std::swap(this->cell_sizes, cell_sizes);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline grid_information<point_t> grid<value_t, point_t, dimensions, rows, columns>::get_grid_information() const
        {
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Grid information not been set."));
            }

            grid_information<point_t> info;
            info.cell_coordinates = this->cell_coordinates;
            info.node_coordinates = this->node_coordinates;
            info.cell_sizes = this->cell_sizes;

            return info;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const typename grid_information<point_t>::array_type& grid<value_t, point_t, dimensions, rows, columns>::get_cell_coordinates() const
        {
#ifdef __tpf_sanity_checks
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell coordinates have not been set."));
            }
#endif

            return this->cell_coordinates;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_cell_coordinates(const coords_t& coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell coordinates have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords[dimension] < 0 || static_cast<std::size_t>(coords[dimension]) < this->extent[dimension].first
                || static_cast<std::size_t>(coords[dimension]) > this->extent[dimension].second)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->cell_coordinates[dimension])[coords[dimension] - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_cell_coordinates(const std::size_t coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell coordinates have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords < 0 || static_cast<std::size_t>(coords) < this->extent[dimension].first
                || static_cast<std::size_t>(coords) > this->extent[dimension].second)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->cell_coordinates[dimension])[coords - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline Eigen::Matrix<point_t, dimensions, 1> grid<value_t, point_t, dimensions, rows, columns>::get_cell_coordinates(const coords_t& coords) const
        {
#ifdef __tpf_sanity_checks
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell coordinates have not been set."));
            }
#endif

            Eigen::Matrix<point_t, dimensions, 1> vector;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                vector[i] = get_cell_coordinates(coords, i);
            }

            return vector;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const typename grid_information<point_t>::array_type& grid<value_t, point_t, dimensions, rows, columns>::get_node_coordinates() const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Node coordinates have not been set."));
            }
#endif

            return this->node_coordinates;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_node_coordinates(const coords_t& coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Node coordinates have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords[dimension] < 0 || static_cast<std::size_t>(coords[dimension]) < this->extent[dimension].first
                || static_cast<std::size_t>(coords[dimension]) > this->extent[dimension].second + 1)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->node_coordinates[dimension])[coords[dimension] - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_node_coordinates(const std::size_t coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Node coordinates have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords < 0 || static_cast<std::size_t>(coords) < this->extent[dimension].first
                || static_cast<std::size_t>(coords) > this->extent[dimension].second + 1)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->node_coordinates[dimension])[coords - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline Eigen::Matrix<point_t, dimensions, 1> grid<value_t, point_t, dimensions, rows, columns>::get_node_coordinates(const coords_t& coords) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Node coordinates have not been set."));
            }
#endif

            Eigen::Matrix<point_t, dimensions, 1> vector;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                vector[i] = get_node_coordinates(coords, i);
            }

            return vector;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const typename grid_information<point_t>::array_type& grid<value_t, point_t, dimensions, rows, columns>::get_cell_sizes() const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell sizes have not been set."));
            }
#endif

            return this->cell_sizes;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_cell_sizes(const coords_t& coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell sizes have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords[dimension] < 0 || static_cast<std::size_t>(coords[dimension]) < this->extent[dimension].first
                || static_cast<std::size_t>(coords[dimension]) > this->extent[dimension].second)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->cell_sizes[dimension])[coords[dimension] - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        const point_t& grid<value_t, point_t, dimensions, rows, columns>::get_cell_sizes(const std::size_t coords, const std::size_t dimension) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell sizes have not been set."));
            }

            if (dimension >= dimensions)
            {
                throw std::runtime_error(__tpf_error_message("Illegal dimension."));
            }
#endif

#ifdef __tpf_range_checks
            if (coords < 0 || static_cast<std::size_t>(coords) < this->extent[dimension].first
                || static_cast<std::size_t>(coords) > this->extent[dimension].second)
            {
                throw std::runtime_error(__tpf_error_message("Access to illegal coordinates."));
            }
#endif

            return (this->cell_sizes[dimension])[coords - this->extent[dimension].first];
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline Eigen::Matrix<point_t, dimensions, 1> grid<value_t, point_t, dimensions, rows, columns>::get_cell_sizes(const coords_t& coords) const
        {
#ifdef __tpf_sanity_checks
            if (!grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Cell sizes have not been set."));
            }
#endif

            Eigen::Matrix<point_t, dimensions, 1> vector;

            for (std::size_t i = 0; i < dimensions; ++i)
            {
                vector[i] = get_cell_sizes(coords, i);
            }

            return vector;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline bool grid<value_t, point_t, dimensions, rows, columns>::has_grid_information() const noexcept
        {
            return this->grid_information_avail;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline constexpr std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_num_components() const noexcept
        {
            return rows * columns;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline constexpr std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_num_rows() const noexcept
        {
            return rows;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline constexpr std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_num_columns() const noexcept
        {
            return columns;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline constexpr std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_dimension() const noexcept
        {
            return dimensions;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_num_elements() const noexcept
        {
            return this->data.get_num_elements();
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline std::size_t grid<value_t, point_t, dimensions, rows, columns>::get_size() const noexcept
        {
            return this->data.get_size();
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const std::vector<value_t>& grid<value_t, point_t, dimensions, rows, columns>::get_data() const noexcept
        {
            return this->data.get_data();
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline std::vector<value_t>& grid<value_t, point_t, dimensions, rows, columns>::get_data() noexcept
        {
            return this->data.get_data();
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline void grid<value_t, point_t, dimensions, rows, columns>::set_data(std::vector<value_t> data_array)
        {
            this->data.set_data(std::move(data_array));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const array<value_t, rows, columns>& grid<value_t, point_t, dimensions, rows, columns>::get_array() const noexcept
        {
            return this->data;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline array<value_t, rows, columns>& grid<value_t, point_t, dimensions, rows, columns>::get_array() noexcept
        {
            return this->data;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline const typename grid<value_t, point_t, dimensions, rows, columns>::return_type grid<value_t, point_t, dimensions, rows, columns>::operator()(const coords_t& coords) const
        {
            return this->data(calculate_index(coords));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline typename grid<value_t, point_t, dimensions, rows, columns>::return_type grid<value_t, point_t, dimensions, rows, columns>::operator()(const coords_t& coords)
        {
            return this->data(calculate_index(coords));
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <typename new_value_t, int new_rows, int new_columns>
        inline grid<new_value_t, point_t, dimensions, new_rows, new_columns> grid<value_t, point_t, dimensions, rows, columns>::copy_structure(std::string name) const
        {
            static_assert(new_rows > 0, "Number of target row components must be larger than zero");
            static_assert(new_columns > 0, "Number of target column components must be larger than zero");

            grid<new_value_t, point_t, dimensions, new_rows, new_columns> copy(std::move(name), this->extent);
            copy.set_grid_information(this->cell_coordinates, this->node_coordinates, this->cell_sizes);

            return copy;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <typename new_value_t, typename new_point_t>
        inline grid<new_value_t, new_point_t, dimensions, rows, columns> grid<value_t, point_t, dimensions, rows, columns>::extract_subgrid(const extent_t& subextent) const
        {
            for (std::size_t d = 0; d < dimensions; ++d)
            {
                if (subextent[d].first < this->extent[d].first)
                {
                    throw std::runtime_error(__tpf_error_message("Subextent must not be larger than the grid."));
                }
            }

            // Split data
            std::size_t num_items = rows * columns;

            for (const auto& ext : subextent)
            {
                num_items *= ext.second - ext.first + 1;
            }

            std::vector<new_value_t> data(num_items);

            for (std::size_t index = 0; index < num_items; ++index)
            {
                // Calculate coords for given index
                Eigen::Matrix<long long, dimensions, 1> coords;

                std::size_t index_rest = index / (rows * columns);
                std::size_t dim = 0;

                for (const auto& ext : subextent)
                {
                    const std::size_t num_items_d = ext.second - ext.first + 1;

                    coords[dim++] = ext.first + index_rest % num_items_d;
                    index_rest /= num_items_d;
                }

                // Copy data
                data[index] = static_cast<new_value_t>(this->data(calculate_index(coords), index % (rows * columns)));
            }

            // Create subgrid
            auto subgrid = grid<new_value_t, new_point_t, dimensions, rows, columns>(get_name(), subextent, std::move(data));

            if (this->grid_information_avail)
            {
                // Split grid information
                auto info = extract_subgrid_information(subextent);

                subgrid.set_grid_information(std::move(info.cell_coordinates), std::move(info.node_coordinates), std::move(info.cell_sizes));
            }

            return subgrid;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <typename new_point_t>
        inline grid_information<new_point_t> grid<value_t, point_t, dimensions, rows, columns>::extract_subgrid_information(const extent_t& subextent) const
        {
            if (!this->grid_information_avail)
            {
                throw std::runtime_error(__tpf_error_message("Grid information has not been set."));
            }

            for (std::size_t d = 0; d < dimensions; ++d)
            {
                if (subextent[d].first < this->extent[d].first)
                {
                    throw std::runtime_error(__tpf_error_message("Subextent must not be larger than the grid."));
                }
            }

            // Split grid information
            grid_information<new_point_t> grid_info;

            grid_info.cell_coordinates.resize(dimensions);
            grid_info.node_coordinates.resize(dimensions);
            grid_info.cell_sizes.resize(dimensions);

            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                grid_info.cell_coordinates[dimension].resize(subextent[dimension].second - subextent[dimension].first + 1);
                grid_info.node_coordinates[dimension].resize(subextent[dimension].second - subextent[dimension].first + 2);
                grid_info.cell_sizes[dimension].resize(subextent[dimension].second - subextent[dimension].first + 1);

                for (std::size_t index = subextent[dimension].first; index <= subextent[dimension].second; ++index)
                {
                    grid_info.cell_coordinates[dimension][index - subextent[dimension].first] = static_cast<new_point_t>(this->cell_coordinates[dimension][index - this->extent[dimension].first]);
                    grid_info.node_coordinates[dimension][index - subextent[dimension].first] = static_cast<new_point_t>(this->node_coordinates[dimension][index - this->extent[dimension].first]);
                    grid_info.cell_sizes[dimension][index - subextent[dimension].first] = static_cast<new_point_t>(this->cell_sizes[dimension][index - this->extent[dimension].first]);
                }

                grid_info.node_coordinates[dimension][subextent[dimension].second - subextent[dimension].first + 1]
                    = static_cast<new_point_t>(this->node_coordinates[dimension][subextent[dimension].second - this->extent[dimension].first + 1]);
            }

            return grid_info;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <typename new_value_t>
        inline void grid<value_t, point_t, dimensions, rows, columns>::fill_subgrid(const grid<new_value_t, point_t, dimensions, rows, columns>& input, const bool fill_grid_information)
        {
            fill_subgrid<new_value_t>(input, input.get_extent(), fill_grid_information);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <typename new_value_t>
        inline void grid<value_t, point_t, dimensions, rows, columns>::fill_subgrid(const grid<new_value_t, point_t, dimensions, rows, columns>& input,
            const extent_t& subextent, const bool fill_grid_information)
        {
            for (std::size_t d = 0; d < dimensions; ++d)
            {
                if (subextent[d].first < input.get_extent()[d].first || subextent[d].first < this->extent[d].first)
                {
                    throw std::runtime_error(__tpf_error_message("Subextent must not be larger than the grid."));
                }
            }

            // Copy data
            std::size_t num_items = rows * columns;

            for (const auto& ext : subextent)
            {
                num_items *= ext.second - ext.first + 1;
            }

            for (std::size_t index = 0; index < num_items; ++index)
            {
                // Calculate coords for given index
                Eigen::Matrix<long long, dimensions, 1> coords;

                std::size_t index_rest = index / (rows * columns);
                std::size_t dim = 0;

                for (const auto& ext : subextent)
                {
                    const std::size_t num_items_d = ext.second - ext.first + 1;

                    coords[dim++] = ext.first + index_rest % num_items_d;
                    index_rest /= num_items_d;
                }

                // Copy data
                this->data(calculate_index(coords), index % (rows * columns)) = static_cast<value_t>(input.get_data()[index]);
            }

            // Copy grid information
            if (fill_grid_information)
            {
                for (std::size_t d = 0; d < dimensions; ++d)
                {
                    for (std::size_t i = subextent[d].first; i <= subextent[d].second; ++i)
                    {
                        this->cell_coordinates[d][i - this->extent[d].first] = input.get_cell_coordinates(i, d);
                        this->node_coordinates[d][i - this->extent[d].first] = input.get_node_coordinates(i, d);
                        this->cell_sizes[d][i - this->extent[d].first] = input.get_cell_sizes(i, d);
                    }

                    this->node_coordinates[d][subextent[d].second + 1 - this->extent[d].first] = input.get_node_coordinates(subextent[d].second + 1, d);
                }
            }
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type
            grid<value_t, point_t, dimensions, rows, columns>::interpolate(const Eigen::Matrix<point_t, dimensions, 1>& point) const
        {
            return interpolate_impl(point);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <bool active>
        inline typename std::enable_if<active, typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type>::type
            grid<value_t, point_t, dimensions, rows, columns>::interpolate_impl(const Eigen::Matrix<point_t, dimensions, 1>& point) const
        {
            using input_type = typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type;

            // Find cell
            const auto cell_opt = find_cell(point);

            if (!cell_opt)
            {
                throw std::runtime_error(__tpf_error_message("Could not find cell for given point."));
            }

            const coords_t cell = *cell_opt;

            // Get all cells, between which the interpolation takes place
            Eigen::Matrix<long long, dimensions, 1> begin, end;

            for (std::size_t d = 0; d < dimensions; ++d)
            {
                coords_t neighbor = cell;

                if (point[d] > this->get_cell_coordinates(cell[d], d))
                {
                    ++neighbor[d];
                }
                else if (point[d] < this->get_cell_coordinates(cell[d], d))
                {
                    --neighbor[d];
                }

                begin[d] = std::min(cell[d], neighbor[d]);
                end[d] = std::max(cell[d], neighbor[d]) + 1;

                if (begin[d] < static_cast<long long>(this->extent[d].first) || end[d] > static_cast<long long>(this->extent[d].second + 1))
                {
                    throw std::runtime_error(__tpf_error_message("Could not get neighborhood for given point. The point is probably located near the boundary of the grid."));
                }
            }

            // Get values and points
            std::vector<std::pair<Eigen::Matrix<point_t, dimensions, 1>, input_type>> point_values;

            algorithm::nested_loop(begin, end, [this](const Eigen::Matrix<long long, dimensions, 1>& coords)
            {
                return std::make_pair(this->get_cell_coordinates(coords), this->operator()(coords));
            }
            , std::back_inserter(point_values));

            // Interpolate
            for (std::size_t d = 0; d < dimensions; ++d)
            {
                // If already finished interpolating, return value
                if (point_values.size() == 1)
                {
                    return point_values[0].second;
                }

                // Interpolate
                std::vector<std::pair<Eigen::Matrix<point_t, dimensions, 1>, input_type>> interpolated_point_values;

                for (std::size_t i = 0; point_values.size() != 0 && i < point_values.size() - 1; i += 2)
                {
                    const auto& left = point_values[i];
                    const auto& right = point_values[i + 1];

                    auto new_position = left.first;
                    new_position[d] = point[d];

                    interpolated_point_values.push_back(std::make_pair(new_position, math::interpolate_linear(point, left.first, right.first, left.second, right.second)));
                }

                point_values = interpolated_point_values;
            }

            // Return interpolated value
            return point_values[0].second;
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        template <bool active>
        inline typename std::enable_if<!active, typename std::conditional<(rows * columns) != 1, Eigen::Matrix<value_t, rows, columns>, value_t>::type>::type
            grid<value_t, point_t, dimensions, rows, columns>::interpolate_impl(const Eigen::Matrix<point_t, dimensions, 1>& point) const
        {
            // Find cell
            const auto cell_opt = find_cell(point);

            if (!cell_opt)
            {
                throw std::runtime_error(__tpf_error_message("Could not find cell for given point."));
            }

            const coords_t cell = *cell_opt;

            // Return value at this cell
            return this->operator()(cell);
        }

        template <typename value_t, typename point_t, int dimensions, int rows, int columns>
        inline std::size_t grid<value_t, point_t, dimensions, rows, columns>::calculate_index(const coords_t& coords) const
        {
            std::size_t index = 0;

            std::size_t multiplier = 1;
            std::size_t coord_index = 0;

            for (auto& ext : this->extent)
            {
                const std::size_t num_elements = ext.second - ext.first + 1;
                const std::size_t coord = coords[coord_index] - ext.first;

#ifdef __tpf_range_checks
                if (coords[coord_index] < 0 || static_cast<std::size_t>(coords[coord_index]) < ext.first || static_cast<std::size_t>(coords[coord_index]) > ext.second)
                {
                    std::stringstream coords_string, extent_string;

                    coords_string << "(";
                    extent_string << "(";

                    for (std::size_t i = 0; i < dimensions - 1; ++i)
                    {
                        coords_string << coords[i] << ", ";
                        extent_string << this->extent[i].first << ":" << this->extent[i].second << ", ";
                    }

                    coords_string << coords[dimensions - 1] << ")";
                    extent_string << this->extent[dimensions - 1].first << ":" << this->extent[dimensions - 1].second << ").\n";

                    throw std::runtime_error(__tpf_error_message("Access to illegal coordinates ", coords_string.str(), " with extent ", extent_string.str()));
                }
#endif

                index += coord * multiplier;
                multiplier *= num_elements;

                ++coord_index;
            }

            return index;
        }
    }
}
