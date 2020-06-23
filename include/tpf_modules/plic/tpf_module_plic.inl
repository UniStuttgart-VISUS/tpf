#include "tpf_module_plic.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_halfspace.h"
#include "tpf/geometry/tpf_intersection.h"
#include "tpf/geometry/tpf_plane.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_polygon.h"
#include "tpf/geometry/tpf_polyhedron.h"

#include "tpf/log/tpf_log.h"

#include "Eigen/Dense"

#include <cmath>
#include <memory>
#include <optional>
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
        inline std::size_t plic<float_t>::get_num_required_ghost_levels()
        {
            return 1;
        }

        template <typename float_t>
        inline plic<float_t>::plic() { }

        template <typename float_t>
        inline std::string plic<float_t>::get_name() const
        {
            return std::string("PLIC");
        }

        template <typename float_t>
        inline void plic<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
            const data::grid<float_t, float_t, 3, 3>& gradients)
        {
            this->fractions = &fractions;
            this->gradients = &gradients;
        }

        template <typename float_t>
        inline void plic<float_t>::set_algorithm_output(data::polydata<float_t>& plic_interface)
        {
            this->plic_interface = &plic_interface;
        }

        template <typename float_t>
        inline void plic<float_t>::set_algorithm_parameters(const std::size_t num_iterations, std::optional<float_t> pertubation)
        {
            this->num_iterations = num_iterations;
            this->perturbation = get_or_default(pertubation, static_cast<float_t>(0.00001L));
        }

        template <typename float_t>
        inline void plic<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;

            data::polydata<float_t>& plic_interface = *this->plic_interface;

            // Calculate PLIC at all interface cells and store the error
#ifdef __tpf_debug
            std::size_t interface_cells = 0;
            std::size_t illegal_cells = 0;

            auto error = std::make_shared<data::array<float_t, 1>>("Error");
#endif

            for (auto z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
            {
                for (auto y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (auto x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                        {
                            // Create PLIC interface
                            auto reconstruction = reconstruct_interface(fractions(coords), gradients(coords),
                                fractions.get_cell_coordinates(coords), fractions.get_cell_sizes(coords), this->num_iterations, this->perturbation);

                            if (reconstruction.first != nullptr)
                            {
                                plic_interface.insert(reconstruction.first);

#ifdef __tpf_debug
                                for (std::size_t i = 0; i < reconstruction.first->get_num_cells(); ++i)
                                {
                                    error->push_back(reconstruction.second);
                                }
#endif
                            }
#ifdef __tpf_debug
                            else
                            {
                                ++illegal_cells;
                            }

                            ++interface_cells;
#endif
                        }
                    }
                }
            }

#ifdef __tpf_debug
            plic_interface.add(error, data::topology_t::CELL_DATA);

            log::info_message(__tpf_info_message("Number of illegal cells: ", illegal_cells, " out of ", interface_cells, " interface cells."));
#endif
        }

        template <typename float_t>
        inline std::pair<std::shared_ptr<geometry::polygon<float_t>>, float_t> plic<float_t>::reconstruct_interface(const float_t vof,
            const Eigen::Matrix<float_t, 3, 1>& gradient, const Eigen::Matrix<float_t, 3, 1>& cell_coordinates,
            const Eigen::Matrix<float_t, 3, 1>& cell_size, const std::size_t num_iterations, float_t perturbation)
        {
            // Auxiliary
            const auto cell_size_half = static_cast<float_t>(0.5) * cell_size;
            const float_t cell_volume = cell_size.prod();

            // Calculate error margin, initial minimum, maximum and ISO value
            const float_t err_margin = static_cast<float_t>(0.0001);

            float_t min_val = err_margin;
            float_t max_val = static_cast<float_t>(1.0) - err_margin;

            // Assume, that this function is only called for cells, where an interface should be calculated.
            // The following definition of starting iso_value will force an interface for the current cell.

            // clamp vof to [min_val, max_val] for numeric stability
            float_t iso_value = std::max(min_val, std::min(max_val, vof));

            // Create cuboid representing the cell
            const auto min_corner = geometry::point<float_t>(cell_coordinates - cell_size_half);
            const auto max_corner = geometry::point<float_t>(cell_coordinates + cell_size_half);

            const geometry::cuboid<float_t> cell(min_corner, max_corner);

            // Get corner point and normal
            const Eigen::Matrix<float_t, 3, 1> corner(
                cell_coordinates[0] + std::copysign(cell_size_half[0], gradient[0]),
                cell_coordinates[1] + std::copysign(cell_size_half[1], gradient[1]),
                cell_coordinates[2] + std::copysign(cell_size_half[2], gradient[2]));

            const Eigen::Matrix<float_t, 3, 1> opposite_corner(
                cell_coordinates[0] - std::copysign(cell_size_half[0], gradient[0]),
                cell_coordinates[1] - std::copysign(cell_size_half[1], gradient[1]),
                cell_coordinates[2] - std::copysign(cell_size_half[2], gradient[2]));

            Eigen::Matrix<float_t, 3, 1> normal = static_cast<float_t>(-1.0) * gradient.normalized();
            // TODO Calculation of orthonormal in plane constructor will throw exception if normal is zero
            //  and __tpf_sanity_checks is defined. Catch this case already here as plane is marked noexcept.
#ifdef __tpf_sanity_checks
            if (normal.isZero())
            {
                throw std::runtime_error(__tpf_error_message("Normal is zero."));
            }
#endif

            // Perform binary search
            std::shared_ptr<geometry::polygon<float_t>> plic = nullptr;
            float_t error = static_cast<float_t>(0.0);

            for (std::size_t i = 0; i < num_iterations; ++i)
            {
                // Calculate PLIC plane for current ISO value
                const geometry::point<float_t> up_point(corner + iso_value * cell_size.norm() * normal);
                const geometry::plane<float_t> plane(up_point, normal);

                // Get intersections and calculate volume
                const std::vector<geometry::point<float_t>> intersections = geometry::template intersect_with<float_t>(plane, cell);

                float_t volume = static_cast<float_t>(0.0);

                if (intersections.size() >= 3)
                {
                    // Add points inside the fluid to the polyhedron and calculate its volume
                    std::vector<geometry::point<float_t>> polyhedron_points(intersections.begin(), intersections.end());

                    for (const auto& corner : cell.get_points())
                    {
                        if (!geometry::template is_in_positive_halfspace<float_t>(geometry::point<float_t>(corner), plane))
                        {
                            polyhedron_points.push_back(geometry::point<float_t>(corner));
                        }
                    }

                    try
                    {
                        volume = geometry::polyhedron<float_t>(polyhedron_points).calculate_volume().get_float_value() / cell_volume;
                    }
                    catch (...)
                    {
                        // Assume co-planar points and thus a minimal velocity
                        volume = static_cast<float_t>(0.0);
                    }
                }
                else
                {
                    // Check where the nearest cell corner is and assign volume accordingly
                    if (geometry::template is_in_positive_halfspace<float_t>(geometry::point<float_t>(corner), plane))
                    {
                        volume = static_cast<float_t>(0.0);
                    }
                    else if (geometry::template is_in_negative_halfspace(geometry::point<float_t>(opposite_corner), plane))
                    {
                        volume = static_cast<float_t>(1.0);
                    }
                    else
                    {
                        log::warning_message(__tpf_warning_message("Degenerate case."));
                        break;
                    }
                }

                // Set new bounds and ISO value
                if (volume > vof)
                {
                    max_val = iso_value;
                }
                else
                {
                    min_val = iso_value;
                }

                iso_value = (max_val + min_val) / 2.0;

                // Create PLIC on last iteration
                if ((error = std::abs(volume - vof)) < err_margin || i == num_iterations - 1)
                {
                    if (intersections.size() >= 3)
                    {
                        // Perturbate points slightly to ensure they lie within the cell cuboid
                        std::vector<geometry::point<float_t>> perturbated_intersections;
                        perturbated_intersections.reserve(intersections.size());

                        perturbation *= cell_size.norm();

                        for (const auto& intersection : intersections)
                        {
                            const math::vec3_t<float_t> perturbated_vertex = intersection.get_vertex() + (cell_coordinates - intersection.get_vertex()).normalized() * perturbation;

                            perturbated_intersections.emplace_back(perturbated_vertex);
                        }

                        plic = std::make_shared<geometry::polygon<float_t>>(perturbated_intersections, true);
                    }

                    break;
                }
            }

            return std::make_pair(plic, error);
        }
    }
}
