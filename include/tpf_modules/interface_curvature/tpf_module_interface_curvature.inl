 #include "tpf_module_interface_curvature.h"

#include "tpf/algorithm/tpf_least_squares.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/math/tpf_vector.h"
#include "tpf/math/tpf_transformer.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t interface_curvature<float_t>::get_num_required_ghost_levels()
        {
            return 1;
        }

        template <typename float_t>
        inline interface_curvature<float_t>::interface_curvature() { }

        template <typename float_t>
        inline std::string interface_curvature<float_t>::get_name() const
        {
            return std::string("Interface Curvature");
        }

        template <typename float_t>
        inline void interface_curvature<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
            const data::grid<float_t, float_t, 3, 3>& gradients, const data::grid<float_t, float_t, 3, 3>& positions)
        {
            this->fractions = &fractions;
            this->gradients = &gradients;
            this->positions = &positions;
        }

        template <typename float_t>
        inline void interface_curvature<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 1>& curvature)
        {
            this->curvature = &curvature;
        }

        template <typename float_t>
        inline void interface_curvature<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;

            data::grid<float_t, float_t, 3, 1>& curvature = *this->curvature;

            // If run on multiple processes, share fractional data
            const data::grid<float_t, float_t, 3, 1> global_fractions = mpi::mpi_grid::all_gather(fractions);

            // Calculate curvature at all cells
            /**
                OpenMP information

                read global:    fractions, gradients, positions, global_fractions
                write global:    curvature
            **/
            #pragma omp parallel for schedule(dynamic) shared(fractions, gradients, positions, curvature)
            for (long long z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                std::size_t z = static_cast<std::size_t>(z_omp);

                for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        curvature(coords) = static_cast<float_t>(0.0L);

                        if (fractions.is_local(coords, get_num_required_ghost_levels()) &&
                            fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                        {
                            const float_t delta = std::min(fractions.get_cell_sizes(coords, 0),
                                std::min(fractions.get_cell_sizes(coords, 1), fractions.get_cell_sizes(coords, 2))) / static_cast<float_t>(2.0L);

                            // Calculate normal and best fitting directions
                            const math::vec3_t<float_t> normal = -gradients(coords).normalized();

                            const auto fitting_directions = math::get_fitting_directions<float_t, 3>(normal);

                            const math::vec3_t<float_t> direction_map(1, 2, 3);
                            const std::array<direction_t, 3> directions = {
                                static_cast<direction_t>(static_cast<int>(fitting_directions[0].dot(direction_map))),
                                static_cast<direction_t>(static_cast<int>(fitting_directions[1].dot(direction_map))),
                                static_cast<direction_t>(static_cast<int>(fitting_directions[2].dot(direction_map))) };

                            // For each direction
                            bool consistent = false;

                            std::vector<math::vec3_t<float_t>> interface_positions;
                            interface_positions.reserve(24);

                            for (std::size_t i = 0; i < 3 && !consistent; ++i)
                            {
                                // Compute the interface height and interface positions
                                const auto curv = height_function_curvature(coords, directions[i], global_fractions);
                                consistent = curv.first;

                                // If the height function returned a valid result, the curvature is set. Else the interface positions returned are added to the set
                                if (!consistent)
                                {
                                    interface_positions.insert(interface_positions.end(), curv.second.second.begin(), curv.second.second.end());
                                }
                                else
                                {
                                    curvature(coords) = curv.second.first;
                                }
                            }

                            if (!consistent)
                            {
                                // If the number of independent positions is smaller than 6
                                std::size_t num_independent_positions = calculate_num_independent_positions(interface_positions, delta);

                                if (num_independent_positions < 6)
                                {
                                    // Replace interface positions with the set of interface positions built from the barycenters of the
                                    // reconstructed interface fragments in a 3 x 3 x 3 stencil
                                    interface_positions.clear();

                                    for (long long i = -1; i <= 1; ++i)
                                    {
                                        for (long long j = -1; j <= 1; ++j)
                                        {
                                            for (long long k = -1; k <= 1; ++k)
                                            {
                                                if (i != 0 || j != 0 || k != 0)
                                                {
                                                    const data::coords3_t neighbor_coords = coords + data::coords3_t(i, j, k);

                                                    if (fractions(neighbor_coords) > static_cast<float_t>(0.0L) && fractions(neighbor_coords) < static_cast<float_t>(1.0L))
                                                    {
                                                        interface_positions.push_back(positions(neighbor_coords));
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    num_independent_positions = calculate_num_independent_positions(interface_positions, delta);
                                }

                                // If the number of independent positions is still too small set a curvature of 0, else apply paraboloid-fitting
                                if (num_independent_positions < 6)
                                {
                                    curvature(coords) = static_cast<float_t>(0.0L);
                                }
                                else
                                {
                                    curvature(coords) = paraboloid_fitted_curvature(normal, interface_positions, positions(coords));
                                }
                            }
                        }
                    }
                }
            }

            mpi::mpi_grid::synchronize_boundaries(curvature);
        }

        template <typename float_t>
        inline float_t interface_curvature<float_t>::calculate_central_derivative(const float_t f0, const float_t f1, const float_t f2, const float_t h0, const float_t h1) const
        {
            const float_t d1 = (f1 - f0) / h0;
            const float_t d2 = (f2 - f1) / h1;

            return (d1 + d2) / static_cast<float_t>(2.0L);
        }

        template <typename float_t>
        inline float_t interface_curvature<float_t>::calculate_central_second_order_derivative(const float_t f0, const float_t f1, const float_t f2, const float_t h0, const float_t h1) const
        {
            const float_t d1 = (f1 - f0) / h0;
            const float_t d2 = (f2 - f1) / h1;

            return (d2 - d1) / ((h0 + h1) / static_cast<float_t>(2.0L));
        }

        template <typename float_t>
        inline std::pair<bool, std::pair<float_t, std::vector<math::vec3_t<float_t>>>> interface_curvature<float_t>::height_function_curvature(const data::coords3_t& coords,
            const direction_t direction, const data::grid<float_t, float_t, 3, 1>& fractions) const
        {
            const std::size_t dim = (direction == direction_t::LEFT || direction == direction_t::RIGHT) ? 0
                : (direction == direction_t::UP || direction == direction_t::DOWN) ? 1 : 2;

            // Calculate interface height and cell coordinates of the base cell
            const auto interface_height_center = calculate_interface_height(coords, direction, fractions);

            // Calculate base cell
            float_t base_center;

            switch (direction) {
            case direction_t::LEFT:
            case direction_t::DOWN:
            case direction_t::BACK:
                base_center = fractions.get_cell_coordinates(interface_height_center.second.first, dim)
                    + fractions.get_cell_sizes(interface_height_center.second.first, dim) / static_cast<float_t>(2.0L);

                break;
            case direction_t::RIGHT:
            case direction_t::UP:
            case direction_t::FRONT:
                base_center = fractions.get_cell_coordinates(interface_height_center.second.first, dim)
                    - fractions.get_cell_sizes(interface_height_center.second.first, dim) / static_cast<float_t>(2.0L);
            }

            // For each of the 8 columns neighboring the central cell
            std::vector<std::pair<bool, std::pair<data::coords3_t, float_t>>> interface_heights;
            interface_heights.reserve(8);

            std::vector<math::vec3_t<float_t>> interface_positions;
            interface_positions.reserve(8);

            data::coords3_t cell_coords;

            bool consistent = true;

            for (long long i = -1; i <= 1; ++i)
            {
                for (long long j = -1; j <= 1; ++j)
                {
                    if (i != 0 || j != 0)
                    {
                        switch (direction) {
                        case direction_t::LEFT:
                        case direction_t::RIGHT:
                            cell_coords[0] = coords[0];
                            cell_coords[1] = coords[1] + i;
                            cell_coords[2] = coords[2] + j;

                            break;
                        case direction_t::UP:
                        case direction_t::DOWN:
                            cell_coords[0] = coords[0] + i;
                            cell_coords[1] = coords[1];
                            cell_coords[2] = coords[2] + j;

                            break;
                        case direction_t::FRONT:
                        case direction_t::BACK:
                        default:
                            cell_coords[0] = coords[0] + i;
                            cell_coords[1] = coords[1] + j;
                            cell_coords[2] = coords[2];
                        }

                        // Calculate interface height and cell coordinates of the base cell
                        auto interface_height = calculate_interface_height(cell_coords, direction, fractions);

                        consistent &= interface_height.first;

                        if (interface_height.first)
                        {
                            math::vec3_t<float_t> interface_position = fractions.get_cell_coordinates(cell_coords);
                            float_t base;

                            // Calculate base cell and interface position
                            switch (direction) {
                            case direction_t::LEFT:
                            case direction_t::DOWN:
                            case direction_t::BACK:
                                base = fractions.get_cell_coordinates(interface_height.second.first, dim)
                                    + fractions.get_cell_sizes(interface_height.second.first, dim) / static_cast<float_t>(2.0L);

                                interface_position[dim] = base - interface_height.second.second;

                                break;
                            case direction_t::RIGHT:
                            case direction_t::UP:
                            case direction_t::FRONT:
                                base = fractions.get_cell_coordinates(interface_height.second.first, dim)
                                    - fractions.get_cell_sizes(interface_height.second.first, dim) / static_cast<float_t>(2.0L);

                                interface_position[dim] = base + interface_height.second.second;
                            }

                            // Set a common origin
                            switch (direction) {
                            case direction_t::LEFT:
                            case direction_t::DOWN:
                            case direction_t::BACK:
                                interface_height.second.second -= base - base_center;

                                break;
                            case direction_t::RIGHT:
                            case direction_t::UP:
                            case direction_t::FRONT:
                            default:
                                interface_height.second.second += base - base_center;
                            }

                            interface_positions.push_back(interface_position);
                        }

                        interface_heights.push_back(interface_height);
                    }
                }
            }

            // If all heights are consistent
            if (consistent)
            {
                // Return the curvature estimated using finite-difference approximations of the derivatives of the discretized height function
                float_t h_x0, h_x1, h_y0, h_y1;

                const data::coords3_t coords_r = coords + data::coords3_t(1, 0, 0);
                const data::coords3_t coords_l = coords + data::coords3_t(-1, 0, 0);
                const data::coords3_t coords_u = coords + data::coords3_t(0, 1, 0);
                const data::coords3_t coords_d = coords + data::coords3_t(0, -1, 0);
                const data::coords3_t coords_f = coords + data::coords3_t(0, 0, 1);
                const data::coords3_t coords_b = coords + data::coords3_t(0, 0, -1);

                switch (direction) {
                case direction_t::LEFT:
                case direction_t::RIGHT:
                    h_x0 = fractions.get_cell_coordinates(coords, 1) - fractions.get_cell_coordinates(coords_d, 1);
                    h_x1 = fractions.get_cell_coordinates(coords_u, 1) - fractions.get_cell_coordinates(coords, 1);
                    h_y0 = fractions.get_cell_coordinates(coords, 2) - fractions.get_cell_coordinates(coords_b, 2);
                    h_y1 = fractions.get_cell_coordinates(coords_f, 2) - fractions.get_cell_coordinates(coords, 2);

                    break;
                case direction_t::UP:
                case direction_t::DOWN:
                    h_x0 = fractions.get_cell_coordinates(coords, 0) - fractions.get_cell_coordinates(coords_l, 0);
                    h_x1 = fractions.get_cell_coordinates(coords_r, 0) - fractions.get_cell_coordinates(coords, 0);
                    h_y0 = fractions.get_cell_coordinates(coords, 2) - fractions.get_cell_coordinates(coords_b, 2);
                    h_y1 = fractions.get_cell_coordinates(coords_f, 2) - fractions.get_cell_coordinates(coords, 2);

                    break;
                case direction_t::FRONT:
                case direction_t::BACK:
                default:
                    h_x0 = fractions.get_cell_coordinates(coords, 0) - fractions.get_cell_coordinates(coords_l, 0);
                    h_x1 = fractions.get_cell_coordinates(coords_r, 0) - fractions.get_cell_coordinates(coords, 0);
                    h_y0 = fractions.get_cell_coordinates(coords, 1) - fractions.get_cell_coordinates(coords_d, 1);
                    h_y1 = fractions.get_cell_coordinates(coords_u, 1) - fractions.get_cell_coordinates(coords, 1);
                }

                // First-order derivatives
                const float_t dx_top = calculate_central_derivative(interface_heights[2].second.second, interface_heights[4].second.second, interface_heights[7].second.second, h_x0, h_x1);
                const float_t dx_center = calculate_central_derivative(interface_heights[1].second.second, interface_height_center.second.second, interface_heights[6].second.second, h_x0, h_x1);
                const float_t dx_bottom = calculate_central_derivative(interface_heights[0].second.second, interface_heights[3].second.second, interface_heights[5].second.second, h_x0, h_x1);

                const float_t dy_center = calculate_central_derivative(interface_heights[3].second.second, interface_height_center.second.second, interface_heights[4].second.second, h_y0, h_y1);

                // Second-order derivatives
                const float_t dxx_center = calculate_central_second_order_derivative(interface_heights[1].second.second, interface_height_center.second.second, interface_heights[6].second.second, h_x0, h_x1);
                const float_t dyy_center = calculate_central_second_order_derivative(interface_heights[3].second.second, interface_height_center.second.second, interface_heights[4].second.second, h_y0, h_y1);

                const float_t dxy_center = calculate_central_derivative(dx_bottom, dx_center, dx_top, h_y0, h_y1);

                const float_t curvature = -((static_cast<float_t>(1.0L) + dx_center * dx_center) * dyy_center
                    - static_cast<float_t>(2.0L) * dx_center * dy_center * dxy_center + (static_cast<float_t>(1.0L) + dy_center * dy_center) * dxx_center)
                    / std::pow(static_cast<float_t>(1.0L) + dx_center * dx_center + dy_center * dy_center, static_cast<float_t>(1.5L));

                return std::make_pair(true, std::make_pair(curvature, interface_positions));
            }
            else
            {
                // Else, return all the interface positions deduced from the consistent heights
                return std::make_pair(false, std::make_pair(static_cast<float_t>(0.0L), interface_positions));
            }
        }

        template <typename float_t>
        inline std::pair<bool, std::pair<data::coords3_t, float_t>> interface_curvature<float_t>::calculate_interface_height(const data::coords3_t& coords,
            const direction_t direction, const data::grid<float_t, float_t, 3, 1>& fractions) const
        {
            const std::size_t dim = (direction == direction_t::LEFT || direction == direction_t::RIGHT) ? 0
                : (direction == direction_t::UP || direction == direction_t::DOWN) ? 1 : 2;

            float_t height = static_cast<float_t>(0.0L);
            data::coords3_t cell_coords;

            // Calculate weighted height of the central cell
            height += fractions(coords) * fractions.get_cell_sizes(coords, dim);

            // Set coordinates to the central cell
            cell_coords = coords;

            // Set interface indicator accordingly
            bool interface = fractions(cell_coords) < static_cast<float_t>(1.0L);

            // While interface = false or at the interface
            while (!interface || (fractions(cell_coords) < static_cast<float_t>(1.0L) && fractions(cell_coords) > static_cast<float_t>(0.0L)))
            {
                // Move along the given direction
                switch (direction)
                {
                case direction_t::LEFT:
                case direction_t::DOWN:
                case direction_t::BACK:
                    --(cell_coords[dim]);

                    break;
                case direction_t::RIGHT:
                case direction_t::UP:
                case direction_t::FRONT:
                    ++(cell_coords[dim]);
                }

                // Check for valid cell coordinates
                if (!fractions.is_on_grid(cell_coords))
                {
                    return std::make_pair(false, std::make_pair(cell_coords, height));
                }

                // Calculate weighted cell height and add to total height
                height += fractions(cell_coords) * fractions.get_cell_sizes(cell_coords, dim);

                // If at the interface, set indicator to true
                if (fractions(cell_coords) < static_cast<float_t>(1.0L) && fractions(cell_coords) > static_cast<float_t>(0.0L))
                {
                    interface = true;
                }
            }

            // If not outside the fluid, return invalid height
            if (fractions(cell_coords) != static_cast<float_t>(0.0L))
            {
                return std::make_pair(false, std::make_pair(cell_coords, height));
            }

            // Set coordinates to the central cell
            cell_coords = coords;

            // Set interface indicator accordingly
            interface = fractions(cell_coords) > static_cast<float_t>(0.0L);

            // While interface = false or at the interface
            while (!interface || (fractions(cell_coords) < static_cast<float_t>(1.0L) && fractions(cell_coords) > static_cast<float_t>(0.0L)))
            {
                // Move contrary to the given direction
                switch (direction)
                {
                case direction_t::LEFT:
                case direction_t::DOWN:
                case direction_t::BACK:
                    ++(cell_coords[dim]);

                    break;
                case direction_t::RIGHT:
                case direction_t::UP:
                case direction_t::FRONT:
                    --(cell_coords[dim]);
                }

                // Check for valid cell coordinates
                if (!fractions.is_on_grid(cell_coords))
                {
                    return std::make_pair(false, std::make_pair(cell_coords, height));
                }

                // Calculate weighted cell height and add to total height
                height += fractions(cell_coords) * fractions.get_cell_sizes(cell_coords, dim);

                // If at the interface, set indicator to true
                if (fractions(cell_coords) < static_cast<float_t>(1.0L) && fractions(cell_coords) > static_cast<float_t>(0.0L))
                {
                    interface = true;
                }
            }

            // If not inside the fluid, return invalid height
            if (fractions(cell_coords) != static_cast<float_t>(1.0L))
            {
                return std::make_pair(false, std::make_pair(cell_coords, height));
            }

            // Return interface height and last cell visited (base cell)
            return std::make_pair(true, std::make_pair(cell_coords, height));
        }

        template <typename float_t>
        inline std::size_t interface_curvature<float_t>::calculate_num_independent_positions(const std::vector<math::vec3_t<float_t>>& positions, const float_t delta) const
        {
            if (positions.size() < 2)
            {
                return positions.size();
            }

            // Create graph with edges between independent positions
            std::vector<std::vector<math::vec3_t<float_t>>> independent_positions_graph(positions.size());

            for (std::size_t i = 0; i < positions.size() - 1; ++i)
            {
                for (std::size_t j = i + 1; j < positions.size(); ++j)
                {
                    if ((positions[i] - positions[j]).norm() >= delta)
                    {
                        independent_positions_graph[i].push_back(positions[j]);
                        independent_positions_graph[j].push_back(positions[i]);
                    }
                }
            }

            // Sort graph for positions with most independent connections
            std::vector<std::size_t> positions_indices(positions.size());
            std::size_t index = 0;

            std::for_each(positions_indices.begin(), positions_indices.end(), [&index](std::size_t& value) { value = index++; });

            std::sort(positions_indices.begin(), positions_indices.end(),
                [&independent_positions_graph](const std::size_t i, const std::size_t j) { return independent_positions_graph[i].size() > independent_positions_graph[j].size(); });

            // Add positions to independent positions list if they have connections to all already inserted elements
            std::vector<math::vec3_t<float_t>> independent_positions;

            independent_positions.push_back(positions[positions_indices[0]]);

            for (std::size_t i = 1; i < positions.size(); ++i)
            {
                bool good = true;

                for (std::size_t j = 0; j < independent_positions.size() && good; ++j)
                {
                    good &= std::find(independent_positions_graph[positions_indices[i]].begin(), independent_positions_graph[positions_indices[i]].end(),
                        independent_positions[j]) != independent_positions_graph[positions_indices[i]].end();
                }

                if (good)
                {
                    independent_positions.push_back(positions[positions_indices[i]]);
                }
            }

            return independent_positions.size();
        }

        template <typename float_t>
        inline float_t interface_curvature<float_t>::paraboloid_fitted_curvature(const math::vec3_t<float_t>& normal,
            const std::vector<math::vec3_t<float_t>>& positions, const math::vec3_t<float_t>& centroid) const
        {
            // Add centroid of the central cell to the interface positions
            std::vector<math::vec3_t<float_t>> positions_transformed;
            positions_transformed.reserve(positions.size() + 1);

            // Transform into coordinate system with centroid as new origin and interface normal as z-axis
            math::transformer<float_t, 3> trafo(centroid, normal);

            for (const auto& position : positions)
            {
                positions_transformed.push_back(trafo.transform_inverse(position));
            }

            positions_transformed.push_back(trafo.transform_inverse(centroid));

            // Fit a paraboloid by minimising F(ai) = sum{1<=j<=n} [zj' - f(ai, xj')]², with f(a_i,x) = a_xx x² + a_yy y² + a_xy xy + a_x x + a_y y + a_1
            const auto poly = algorithm::least_squares(positions_transformed.begin(), positions_transformed.end());

            // Return the mean curvature at the origin: K = -2 * [a_xx (1 + a_y²) + a_yy (1 + a_x²) - a_xy a_x a_y] / [(1 + a_x² + a_y²)^(3/2)]
            const float_t numerator = poly.a_xx * (static_cast<float_t>(1.0L) + poly.a_y * poly.a_y) + poly.a_yy
                * (static_cast<float_t>(1.0L) + poly.a_x * poly.a_x) - poly.a_xy * poly.a_x * poly.a_y;
            const float_t denominator = std::pow(static_cast<float_t>(1.0L) + poly.a_x * poly.a_x + poly.a_y * poly.a_y, static_cast<float_t>(1.5L));

            return -static_cast<float_t>(2.0L) * numerator / denominator;
        }
    }
}
