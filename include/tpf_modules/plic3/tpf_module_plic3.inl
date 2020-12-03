#include "tpf_module_plic3.h"

#include "../interface_gradient/tpf_module_interface_gradient.h"
#include "../plic/tpf_module_plic.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"
#include "tpf/data/tpf_stencil.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_distance.h"
#include "tpf/geometry/tpf_halfspace.h"
#include "tpf/geometry/tpf_intersection.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_vector.h"

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
        inline plic3<float_t>::plic3() { }

        template <typename float_t>
        inline std::string plic3<float_t>::get_name() const
        {
            return std::string("PLIC3");
        }

        template <typename float_t>
        inline void plic3<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& f,
            const data::grid<float_t, float_t, 3, 1>& f3, const data::grid<float_t, float_t, 3, 3>& f_norm_3ph)
        {
            this->f = &f;
            this->f3 = &f3;
            this->f_norm_3ph = &f_norm_3ph;
        }

        template <typename float_t>
        inline void plic3<float_t>::set_algorithm_output(data::polydata<float_t>& plic3_interface)
        {
            this->plic3_interface = &plic3_interface;
        }

        template <typename float_t>
        inline void plic3<float_t>::set_algorithm_parameters(float_t epsilon, float_t error_margin_plic,
            float_t error_margin_plic3, std::size_t num_iterations_plic, std::size_t num_iterations_plic3,
            std::optional<float_t> perturbation, std::optional<bool> hide_error_cubes)
        {
            this->epsilon_ = epsilon;
            this->error_margin_plic_ = error_margin_plic;
            this->error_margin_plic3_ = error_margin_plic3;
            this->num_iterations_plic_ = num_iterations_plic;
            this->num_iterations_plic3_ = num_iterations_plic3;
            this->perturbation_ = get_or_default<float_t>(perturbation, static_cast<float_t>(0.00001L));
            this->hide_error_cubes_ = get_or_default<float_t>(hide_error_cubes, false);
        }

        template <typename float_t>
        inline void plic3<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& f = *this->f;
            const data::grid<float_t, float_t, 3, 1>& f3 = *this->f3;
            const data::grid<float_t, float_t, 3, 3>& f_norm_3ph = *this->f_norm_3ph;

            data::polydata<float_t>& plic3_interface = *this->plic3_interface;
            auto types = std::make_shared<data::array<int, 1>>("type");
            auto gradients_f = std::make_shared<data::array<float_t, 3>>("grad");
            auto gradients_f3 = std::make_shared<data::array<float_t, 3>>("grad3");
            auto values_f = std::make_shared<data::array<float_t, 1>>("f");
            auto values_f3 = std::make_shared<data::array<float_t, 1>>("f3");
            auto values_n3 = std::make_shared<data::array<float_t, 3>>("fnorm3ph");

            const float_t epsilon = this->epsilon_;
            const float_t epsilon_one = static_cast<float_t>(1.0L) - epsilon;

            const auto extent = f.get_extent();

            for (auto z = extent[2].first; z <= extent[2].second; ++z) {
                for (auto y = extent[1].first; y <= extent[1].second; ++y) {
                    for (auto x = extent[0].first; x <= extent[0].second; ++x) {
                        const data::coords3_t coords(x, y, z);

                        auto f_value = f(coords);
                        auto f3_value = f3(coords);
                        Eigen::Matrix<float_t, 3, 1> n3_value = f_norm_3ph(coords);
                        Eigen::Matrix<float_t, 3, 1> n3_output = Eigen::Matrix<float_t, 3, 1>::Zero();

                        // empty cells
                        if (f_value < epsilon && f3_value < epsilon) {
                            continue;
                        }
                        // full solid cell - no interface
                        if (f3_value > epsilon_one && !this->is_interface(f3, coords, epsilon))
                        {
                            continue;
                        }
                        // full fluid cell - no interface
                        if (f_value > epsilon_one && !this->is_interface(f, coords, epsilon))
                        {
                            continue;
                        }

                        bool done = false;
                        polygon_t type = polygon_t::ERROR;
                        math::vec3_t<float_t> grad = interface_gradient<float_t>::calculate_gradient(coords, f); // TODO perfomance, not needed in all cases below
                        math::vec3_t<float_t> grad3 = interface_gradient<float_t>::calculate_gradient(coords, f3); // TODO perfomance, not needed in all cases below

                        // full solid cell - interface
                        if (!done && f3_value > epsilon_one)
                        {
                            auto reconstruction = plic<float_t>::reconstruct_interface(epsilon_one, grad3, f3.get_cell_coordinates(coords), f3.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                auto neighbor = neighbor_type(f3, f, coords, epsilon);
                                switch (neighbor) {
                                    case neighbor_t::MIXED: {
                                        type = polygon_t::SOLID_MIX;
                                        break;
                                    }
                                    case neighbor_t::EMPTY: {
                                        type = polygon_t::SOLID_GAS;
                                        break;
                                    }
                                    case neighbor_t::OTHER: {
                                        type = polygon_t::SOLID_FLUID;
                                        break;
                                    }
                                    default: {
                                        type = polygon_t::ERROR;
                                        break;
                                    }
                                }

                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // full fluid cell - interface
                        if (!done && f_value > epsilon_one)
                        {
                            auto neighbor = neighbor_type(f, f3, coords, epsilon);
                            if (neighbor == neighbor_t::OTHER) {
                                continue;
                            }
                            auto reconstruction = plic<float_t>::reconstruct_interface(epsilon_one, grad, f.get_cell_coordinates(coords), f.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                type = polygon_t::FLUID_GAS;
                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // only solid interface
                        if (!done && f_value < epsilon)
                        {
                            auto reconstruction = plic<float_t>::reconstruct_interface(f3_value, grad3, f3.get_cell_coordinates(coords), f3.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                type = polygon_t::SOLID_GAS;
                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // only fluid interface - without solid neighbor
                        if (!done && f3_value < epsilon && !this->has_non_zero_neighbor(f3, coords, epsilon))
                        {
                            auto reconstruction = plic<float_t>::reconstruct_interface(f_value, grad, f.get_cell_coordinates(coords), f.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                type = polygon_t::FLUID_GAS;
                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // only solid + fluid
                        if (!done && f3_value + f_value > epsilon_one)
                        {
                            auto reconstruction = plic<float_t>::reconstruct_interface(f3_value, grad3, f3.get_cell_coordinates(coords), f3.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                type = polygon_t::SOLID_FLUID;
                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // only fluid interface, but neighbor cell with solid (Three-Phase - Typ 2)
                        if (!done && f3_value < epsilon)
                        {
                            if (std::isnan(n3_value[0]) || std::isnan(n3_value[1]) || std::isnan(n3_value[2]))
                            {
                                throw std::runtime_error(__tpf_error_message("Normal contains nan values."));
                            }
                            auto reconstruction = plic<float_t>::reconstruct_interface(f_value, -n3_value, f.get_cell_coordinates(coords), f.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (reconstruction.first != nullptr)
                            {
                                type = polygon_t::FLUID_TYPE2;
                                n3_output = n3_value;
                                plic3_interface.insert(reconstruction.first);
                                done = true;
                            }
                        }
                        // three pase cell (Three-Phase - Typ 1 or 3)
                        if (!done && f3_value >= epsilon && f_value >= epsilon && f3_value + f_value <= epsilon_one)
                        {
                            if (std::isnan(n3_value[0]) || std::isnan(n3_value[1]) || std::isnan(n3_value[2]))
                            {
                                throw std::runtime_error(__tpf_error_message("Normal contains nan values."));
                            }
                            auto f3_reconstruction = plic<float_t>::reconstruct_interface(f3_value, grad3, f3.get_cell_coordinates(coords), f3.get_cell_sizes(coords), this->error_margin_plic_, this->num_iterations_plic_, this->perturbation_);
                            if (f3_reconstruction.first != nullptr)
                            {
                                // Add extra polygon, therefore directly add to lists and not set done=true
                                plic3_interface.insert(f3_reconstruction.first);
                                types->push_back(polygon_t::SOLID_MIX);
                                gradients_f->push_back(grad);
                                gradients_f3->push_back(grad3);
                                values_f->push_back(f_value);
                                values_f3->push_back(f3_value);
                                values_n3->push_back(n3_output); // Still want zero n3_value here.

                                // interface plane
                                const auto& f3_interface = f3_reconstruction.first;
                                const geometry::plane<float_t> f3_plane(f3_interface->get_points()[0], grad3);

                                // the cell
                                const auto& cell_coordinates = f3.get_cell_coordinates(coords);
                                const auto& cell_size = f3.get_cell_sizes(coords);
                                const auto cell_size_half = static_cast<float_t>(0.5) * cell_size;

                                // Create cuboid representing the cell
                                const auto min_corner = geometry::point<float_t>(cell_coordinates - cell_size_half);
                                const auto max_corner = geometry::point<float_t>(cell_coordinates + cell_size_half);
                                const geometry::cuboid<float_t> cell(min_corner, max_corner);

                                // Generate list of polyhedron points
                                // add interface points to polyhedron
                                const auto& f3_interface_points = f3_interface->get_points();
                                std::vector<geometry::point<float_t>> polyhedron_points(f3_interface_points.begin(), f3_interface_points.end());
                                // add cell points to polyhedron
                                for (const auto& corner : cell.get_points())
                                {
                                    if (!geometry::template is_in_positive_halfspace<float_t>(geometry::point<float_t>(corner), f3_plane))
                                    {
                                        polyhedron_points.push_back(geometry::point<float_t>(corner));
                                    }
                                }
                                geometry::polyhedron<float_t> poly_cell(polyhedron_points);

                                auto f_reconstruction = reconstruct_f3_interface(poly_cell, f_value, n3_value, f3.get_cell_coordinates(coords), f3.get_cell_sizes(coords), this->error_margin_plic3_, this->num_iterations_plic3_, this->perturbation_);
                                if (f_reconstruction.first != nullptr)
                                {
                                    type = polygon_t::FLUID_TYPE13;
                                    n3_output = n3_value;
                                    plic3_interface.insert(f_reconstruction.first);
                                    done = true;
                                }
                            }
                        }

                        if (!done && !this->hide_error_cubes_) {
                            // draw error cube
                            const auto pos = f.get_cell_coordinates(coords);
                            const auto size = f.get_cell_sizes(coords);

                            auto p_min = tpf::geometry::point<float_t>(pos - 0.25 * size);
                            auto p_max = tpf::geometry::point<float_t>(pos + 0.25 * size);
                            auto cube = std::make_shared<tpf::geometry::cuboid<float_t>>(p_min, p_max);

                            plic3_interface.insert(cube);
                            done = true;
                        }
                        if (done) {
                            types->push_back(type);
                            gradients_f->push_back(grad);
                            gradients_f3->push_back(grad3);
                            values_f->push_back(f_value);
                            values_f3->push_back(f3_value);
                            values_n3->push_back(n3_output);
                        }
                    }
                }
            }

            plic3_interface.add(types, tpf::data::topology_t::OBJECT_DATA);
            plic3_interface.add(gradients_f, tpf::data::topology_t::OBJECT_DATA);
            plic3_interface.add(gradients_f3, tpf::data::topology_t::OBJECT_DATA);
            plic3_interface.add(values_f, tpf::data::topology_t::OBJECT_DATA);
            plic3_interface.add(values_f3, tpf::data::topology_t::OBJECT_DATA);
            plic3_interface.add(values_n3, tpf::data::topology_t::OBJECT_DATA);
        }

        template <typename float_t>
        inline bool plic3<float_t>::is_interface(const data::grid<float_t, float_t, 3, 1>& f, const data::coords3_t& coords, float_t epsilon) const
        {
            using stencil_t = data::stencil<const float_t, float_t, 3, 1>;
            const stencil_t stencil(f, coords, 3, stencil_t::behavior_t::REPEAT);

            const data::coords3_t neighbors[6] = {
                data::coords3_t(0, 1, 1),
                data::coords3_t(2, 1, 1),
                data::coords3_t(1, 0, 1),
                data::coords3_t(1, 2, 1),
                data::coords3_t(1, 1, 0),
                data::coords3_t(1, 1, 2)
            };

            for (const auto& neighbor : neighbors)
            {
                if (stencil(neighbor) < epsilon)
                {
                    return true;
                }
            }
            return false;
        }

        template <typename float_t>
        inline typename plic3<float_t>::neighbor_t plic3<float_t>::neighbor_type(const data::grid<float_t, float_t, 3, 1>& f1,
                                                                                 const data::grid<float_t, float_t, 3, 1>& f2,
                                                                                 const data::coords3_t& coords,
                                                                                 float_t epsilon) const
        {
            using stencil_t = data::stencil<const float_t, float_t, 3, 1>;
            const stencil_t stencil1(f1, coords, 3, stencil_t::behavior_t::REPEAT);
            const stencil_t stencil2(f2, coords, 3, stencil_t::behavior_t::REPEAT);

            bool other_neighbor = false;
            bool empty_neighbor = false;

            const data::coords3_t neighbors[6] = {
                data::coords3_t(0, 1, 1),
                data::coords3_t(2, 1, 1),
                data::coords3_t(1, 0, 1),
                data::coords3_t(1, 2, 1),
                data::coords3_t(1, 1, 0),
                data::coords3_t(1, 1, 2)
            };

            for (const auto& neighbor : neighbors) {
                auto s1 = stencil1(neighbor);
                auto s2 = stencil2(neighbor);
                if (s2 > epsilon) {
                    other_neighbor = true;
                }
                if (s1 + s2 < static_cast<float_t>(1.0) - epsilon) {
                    empty_neighbor = true;
                }
            }

            if (other_neighbor && empty_neighbor) {
                return neighbor_t::MIXED;
            }
            if (other_neighbor) {
                return neighbor_t::OTHER;
            }
            if (empty_neighbor) {
                return neighbor_t::EMPTY;
            }
            return neighbor_t::SAME;
        }

        template <typename float_t>
        inline bool plic3<float_t>::has_non_zero_neighbor(const data::grid<float_t, float_t, 3, 1>& f, const data::coords3_t& coords, float_t epsilon) const
        {
            using stencil_t = data::stencil<const float_t, float_t, 3, 1>;
            const stencil_t stencil(f, coords, 3, stencil_t::behavior_t::REPEAT);

            for (std::size_t i = 0; i < 2; ++i) {
                for (std::size_t j = 0; j < 2; ++j) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        if (i == 1 && j == 1 && k == 1) {
                            continue;
                        }
                        if (stencil(data::coords3_t(i, j, k)) > epsilon) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        template <typename float_t>
        inline std::pair<std::shared_ptr<geometry::polygon<float_t>>, float_t> plic3<float_t>::reconstruct_f3_interface(const geometry::polyhedron<float_t>& poly_cell,
            float_t vof, const Eigen::Matrix<float_t, 3, 1>& norm, const Eigen::Matrix<float_t, 3, 1>& cell_coordinates,
            const Eigen::Matrix<float_t, 3, 1>& cell_size, float_t error_margin, std::size_t num_iterations, float_t perturbation)
        {
            // Calculate error margin, initial minimum, maximum and ISO value
            float_t min_val = error_margin;
            float_t max_val = static_cast<float_t>(1.0) - error_margin;

            if (vof <= min_val || vof >= max_val)
            {
                return std::make_pair(nullptr, error_margin);
            }

            float_t iso_value = vof;

            // points
            const auto& points = poly_cell.get_points();
            const geometry::plane<float_t> dist_plane(points[0], norm);

            // determine corner and opposite corner
            const auto& poly_points = poly_cell.get_points();
            float_t min_dist = 0;
            float_t max_dist = 0;
            size_t liquid_corner_idx = 0;
            size_t vaporous_corner_idx = 0;

            for (std::size_t i = 0; i < poly_points.size(); ++i) {
                const auto point = geometry::point<float_t>(poly_points[i]);
                float_t factor = geometry::is_in_positive_halfspace(point, dist_plane) ? static_cast<float_t>(1.0) : static_cast<float_t>(-1.0);
                float_t dist = factor * geometry::calculate_distance<float_t>(point, dist_plane);
                if (dist < min_dist) {
                    min_dist = dist;
                    liquid_corner_idx = i;
                }
                if (dist > max_dist) {
                    max_dist = dist;
                    vaporous_corner_idx = i;
                }
            }

            // Get corner point and normal
            const Eigen::Matrix<float_t, 3, 1> corner = poly_points[liquid_corner_idx];

            const Eigen::Matrix<float_t, 3, 1> opposite_corner = poly_points[vaporous_corner_idx];

            // Perform binary search
            std::shared_ptr<geometry::polygon<float_t>> plic = nullptr;
            float_t error = static_cast<float_t>(0.0);

            const auto poly_cell_diag = (opposite_corner - corner).norm();
            const auto poly_cell_volume = poly_cell.calculate_volume();
            const auto total_cell_diag = cell_size.norm();
            const auto total_cell_volume = cell_size.prod();

            for (std::size_t i = 0; i < num_iterations; ++i)
            {
                // Calculate PLIC plane for current ISO value
                const geometry::point<float_t> up_point(corner + iso_value * poly_cell_diag * norm);
                const geometry::plane<float_t> plane(up_point, norm);

                // Get intersections and calculate volume
                const std::vector<geometry::point<float_t>> intersections = geometry::template intersect_with<float_t>(plane, poly_cell);

                float_t volume = static_cast<float_t>(0.0);

                if (intersections.size() >= 3)
                {
                    // Add points inside the fluid to the polyhedron and calculate its volume
                    std::vector<geometry::point<float_t>> polyhedron_points(intersections.begin(), intersections.end());

                    for (const auto& corner : poly_cell.get_points())
                    {
                        if (!geometry::template is_in_positive_halfspace<float_t>(geometry::point<float_t>(corner), plane))
                        {
                            polyhedron_points.push_back(geometry::point<float_t>(corner));
                        }
                    }

                    try
                    {
                        volume = geometry::polyhedron<float_t>(polyhedron_points).calculate_volume().get_float_value() / total_cell_volume;
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
                if ((error = std::abs(volume - vof)) < error_margin || i == num_iterations - 1)
                {
                    if (intersections.size() >= 3)
                    {
                        // Perturbate points slightly to ensure they lie within the cell cuboid
                        std::vector<geometry::point<float_t>> perturbated_intersections;
                        perturbated_intersections.reserve(intersections.size());

                        perturbation *= poly_cell_diag;

                        for (const auto& intersection : intersections)
                        {
                            const math::vec3_t<float_t> perturbated_vertex = intersection.get_vertex() + (cell_coordinates - intersection.get_vertex()).normalized() * perturbation;

                            perturbated_intersections.emplace_back(perturbated_vertex);
                        }

                        plic = std::make_shared<geometry::polygon<float_t>>(perturbated_intersections, norm, true);
                    }

                    break;
                }
            }

            return std::make_pair(plic, error);
        }
    }
}
