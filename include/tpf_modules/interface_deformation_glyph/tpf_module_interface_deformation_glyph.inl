#include "tpf_module_interface_deformation_glyph.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"
#include "tpf/data/tpf_polydata_element.h"

#include "tpf/geometry/tpf_mesh.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_math_defines.h"
#include "tpf/math/tpf_transformer.h"
#include "tpf/math/tpf_vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t interface_deformation_glyph<float_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t>
        inline interface_deformation_glyph<float_t>::interface_deformation_glyph() { }

        template <typename float_t>
        inline std::string interface_deformation_glyph<float_t>::get_name() const
        {
            return std::string("Interface Deformation Glyph");
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::set_algorithm_input(
            const data::grid<float_t, float_t, 3, 1>& vof,
            const data::grid<float_t, float_t, 3, 3>& positions,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> gradients,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> velocities,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_min,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_max,
            opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_min,
            opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_max,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_min,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_max,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_polynomial)
        {
            this->vof = &vof;
            this->positions = &positions;

            this->gradients = get_or_default<const data::grid<float_t, float_t, 3, 3>>(gradients);

            this->velocities = get_or_default<const data::grid<float_t, float_t, 3, 3>>(velocities);

            this->stretching_min = get_or_default<const data::grid<float_t, float_t, 3, 3>>(stretching_min);
            this->stretching_max = get_or_default<const data::grid<float_t, float_t, 3, 3>>(stretching_max);

            this->bending_min = get_or_default<const data::grid<float_t, float_t, 3, 1>>(bending_min);
            this->bending_max = get_or_default<const data::grid<float_t, float_t, 3, 1>>(bending_max);
            this->bending_direction_min = get_or_default<const data::grid<float_t, float_t, 3, 3>>(bending_direction_min);
            this->bending_direction_max = get_or_default<const data::grid<float_t, float_t, 3, 3>>(bending_direction_max);
            this->bending_polynomial = get_or_default<const data::grid<float_t, float_t, 3, 3>>(bending_polynomial);
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::set_algorithm_output(data::polydata<float_t>& velocity_glyphs,
            std::array<std::reference_wrapper<data::polydata<float_t>>, 3> stretching_glyphs,
            std::array<std::reference_wrapper<data::polydata<float_t>>, 2> bending_glyphs)
        {
            this->velocity_glyphs = &velocity_glyphs;

            this->stretching_glyph_rings = &stretching_glyphs[0].get();
            this->stretching_glyph_references = &stretching_glyphs[1].get();
            this->stretching_glyph_strips = &stretching_glyphs[2].get();

            this->bending_glyph_discs = &bending_glyphs[0].get();
            this->bending_glyph_strips = &bending_glyphs[1].get();
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::set_algorithm_parameters(const bool velocity_glyph,
            const bool stretching_glyph, const bool bending_glyph, const float_t timestep,
            interface_deformation_glyph_aux::velocity_params_t<float_t> velocity_parameters,
            interface_deformation_glyph_aux::stretching_params_t<float_t> stretching_parameters,
            interface_deformation_glyph_aux::bending_params_t<float_t> bending_parameters)
        {
            this->velocity_glyph = velocity_glyph;
            this->stretching_glyph = stretching_glyph;
            this->bending_glyph = bending_glyph;

            this->timestep = timestep;

            this->velocity_parameters = velocity_parameters;
            this->stretching_parameters = stretching_parameters;
            this->bending_parameters = bending_parameters;
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;

            // Calculate average cell size
            const auto size_x = (vof.get_node_coordinates()[0].back() - vof.get_node_coordinates()[0].front())
                / (vof.get_extent()[0].second - vof.get_extent()[0].first + 1);
            const auto size_y = (vof.get_node_coordinates()[1].back() - vof.get_node_coordinates()[1].front())
                / (vof.get_extent()[1].second - vof.get_extent()[1].first + 1);
            const auto size_z = (vof.get_node_coordinates()[2].back() - vof.get_node_coordinates()[2].front())
                / (vof.get_extent()[2].second - vof.get_extent()[2].first + 1);

            const auto average_cell_size = (size_x + size_y + size_z) / static_cast<float_t>(3.0);

            // Create velocity glyph
            if (this->velocity_glyph && this->velocities != nullptr)
            {
                instantiate_velocity_glyphs(create_velocity_glyph_template(
                    this->velocity_parameters.resolution, this->velocity_parameters.shaft_length, this->velocity_parameters.shaft_thickness),
                    average_cell_size, this->velocity_parameters.size, this->velocity_parameters.scalar, this->velocity_parameters.fixed_scalar);
            }

            // Create stretching glyph
            if (this->stretching_glyph && this->gradients != nullptr && this->stretching_min != nullptr && this->stretching_max != nullptr)
            {
                instantiate_stretching_glpyh(create_stretching_glyph_template(this->stretching_parameters.disc_resolution,
                    this->stretching_parameters.disc_bending, this->stretching_parameters.hole_radius,
                    this->stretching_parameters.offset, this->stretching_parameters.strip_size, this->stretching_parameters.reference_size,
                    this->stretching_parameters.z_offset), average_cell_size, this->stretching_parameters.size_scalar,
                    this->stretching_parameters.scalar, this->stretching_parameters.sharpness,
                    this->stretching_parameters.show_strip, this->stretching_parameters.show_reference);
            }

            // Create bending glyph
            if (this->bending_glyph && this->gradients != nullptr && this->bending_min != nullptr && this->bending_max != nullptr &&
                this->bending_direction_min != nullptr && this->bending_direction_max != nullptr)
            {
                instantiate_bending_glyph(create_bending_glyph_template(this->bending_parameters.disc_resolution,
                    this->bending_parameters.polynomial_resolution, this->bending_parameters.offset,
                    this->bending_parameters.strip_size, this->bending_parameters.z_offset), average_cell_size,
                    this->bending_parameters.size_scalar, this->bending_parameters.scalar, this->bending_parameters.show_strip);
            }
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_velocity_glyph_template(std::size_t resolution,
            float_t shaft_length, float_t shaft_thickness) const -> velocity_glyph_t
        {
            // Clamp parameter values
            resolution = std::max(resolution, static_cast<std::size_t>(3));
            shaft_thickness = std::clamp(shaft_thickness, static_cast<float_t>(0.01), static_cast<float_t>(0.5));
            shaft_length = std::clamp(shaft_length, static_cast<float_t>(0.02), static_cast<float_t>(0.98));

            const auto tip_thickness = static_cast<float_t>(2.0) * shaft_thickness;

            const auto increment = static_cast<float_t>((2.0 * math::pi<float_t>) / resolution);
            const auto offset = 0.01;

            // Create arrow glyph...
            velocity_glyph_t arrow_glyph = std::make_shared<geometry::mesh<float_t>>();

            {
                // ... shaft
                {
                    std::vector<CGAL::SM_Vertex_index> origin_1_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> origin_2_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> target_indices(resolution);

                    const auto cap_center = arrow_glyph->add_point(geometry::point<float_t>(0.0, 0.0, 0.0));

                    origin_1_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        0.0,
                        shaft_thickness * std::cos(static_cast<float_t>(0.0)),
                        shaft_thickness * std::sin(static_cast<float_t>(0.0))));

                    origin_2_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        0.0,
                        shaft_thickness * std::cos(static_cast<float_t>(0.0)),
                        shaft_thickness * std::sin(static_cast<float_t>(0.0))));

                    target_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        shaft_length,
                        shaft_thickness * std::cos(static_cast<float_t>(0.0)),
                        shaft_thickness * std::sin(static_cast<float_t>(0.0))));

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        origin_1_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            0.0,
                            shaft_thickness * std::cos(angle),
                            shaft_thickness * std::sin(angle)));

                        origin_2_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            0.0,
                            shaft_thickness * std::cos(angle),
                            shaft_thickness * std::sin(angle)));

                        target_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            shaft_length,
                            shaft_thickness * std::cos(angle),
                            shaft_thickness * std::sin(angle)));

                        arrow_glyph->add_face(cap_center, origin_1_indices[i - 1], origin_1_indices[i]);
                        arrow_glyph->add_face({ origin_2_indices[i - 1], target_indices[i - 1], target_indices[i], origin_2_indices[i] });
                    }

                    arrow_glyph->add_face(cap_center, origin_1_indices.back(), origin_1_indices.front());
                    arrow_glyph->add_face({ origin_2_indices.back(), target_indices.back(), target_indices.front(), origin_2_indices.front() });
                }

                // ... tip
                {
                    std::vector<CGAL::SM_Vertex_index> origin_1_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> origin_2_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> target_1_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> target_2_indices(resolution);

                    const auto tip_center = arrow_glyph->add_point(geometry::point<float_t>(1.0, 0.0, 0.0));
                    const auto tip_cap_center = arrow_glyph->add_point(geometry::point<float_t>(shaft_length, 0.0, 0.0));

                    origin_1_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        shaft_length,
                        tip_thickness * std::cos(static_cast<float_t>(0.0)),
                        tip_thickness * std::sin(static_cast<float_t>(0.0))));

                    origin_2_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        shaft_length,
                        tip_thickness * std::cos(static_cast<float_t>(0.0)),
                        tip_thickness * std::sin(static_cast<float_t>(0.0))));

                    target_1_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        1.0 - offset,
                        tip_thickness * (offset / (1.0 - shaft_length)) * std::cos(static_cast<float_t>(0.0)),
                        tip_thickness * (offset / (1.0 - shaft_length)) * std::sin(static_cast<float_t>(0.0))));

                    target_2_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        1.0 - offset,
                        tip_thickness * (offset / (1.0 - shaft_length)) * std::cos(static_cast<float_t>(0.0)),
                        tip_thickness * (offset / (1.0 - shaft_length)) * std::sin(static_cast<float_t>(0.0))));

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        origin_1_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            shaft_length,
                            tip_thickness * std::cos(angle),
                            tip_thickness * std::sin(angle)));

                        origin_2_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            shaft_length,
                            tip_thickness * std::cos(angle),
                            tip_thickness * std::sin(angle)));

                        target_1_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            1.0 - offset,
                            tip_thickness * (offset / (1.0 - shaft_length)) * std::cos(angle),
                            tip_thickness * (offset / (1.0 - shaft_length)) * std::sin(angle)));

                        target_2_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            1.0 - offset,
                            tip_thickness * (offset / (1.0 - shaft_length)) * std::cos(angle),
                            tip_thickness * (offset / (1.0 - shaft_length)) * std::sin(angle)));

                        arrow_glyph->add_face(tip_cap_center, origin_1_indices[i - 1], origin_1_indices[i]);
                        arrow_glyph->add_face({ origin_2_indices[i - 1], target_1_indices[i - 1], target_1_indices[i], origin_2_indices[i] });
                        arrow_glyph->add_face(tip_center, target_2_indices[i], target_2_indices[i - 1]);
                    }

                    arrow_glyph->add_face(tip_cap_center, origin_1_indices.back(), origin_1_indices.front());
                    arrow_glyph->add_face({ origin_2_indices.back(), target_1_indices.back(), target_1_indices.front(), origin_2_indices.front() });
                    arrow_glyph->add_face(tip_center, target_2_indices.front(), target_2_indices.back());
                }
            }

            // Create second arrow glyph
            const math::transformer<float, 3> trafo(
                Eigen::Matrix<float_t, 3, 1>(-1.0, 0.0, 0.0),
                Eigen::Matrix<float_t, 3, 1>(1.0, 0.0, 0.0),
                Eigen::Matrix<float_t, 3, 1>(0.0, 1.0, 0.0),
                Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0));

            return arrow_glyph->merge(*std::static_pointer_cast<geometry::mesh<float_t>>(arrow_glyph->clone(trafo)));
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::instantiate_velocity_glyphs(const velocity_glyph_t& glyph_template,
            const float_t average_cell_size, const interface_deformation_glyph_aux::arrow_size_t arrow_size,
            float_t arrow_scalar, float_t arrow_fixed_scalar)
        {
            // Clamp parameter values
            arrow_scalar = std::max(arrow_scalar, static_cast<float_t>(0.00001 * average_cell_size));
            arrow_fixed_scalar = std::max(arrow_fixed_scalar, static_cast<float_t>(0.00001 * average_cell_size));

            // Get input data
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& velocities = *this->velocities;

            // Create output data arrays
            auto magnitudes = std::make_shared<data::array<float_t>>("Velocity Magnitude");

            // Iterate over all interface cells
            #pragma omp parallel for schedule(dynamic) default(none) shared(arrow_scalar, arrow_fixed_scalar)
            for (long long z = static_cast<long long>(vof.get_extent()[2].first);
                z <= static_cast<long long>(vof.get_extent()[2].second); ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (vof(coords) > 0 && vof(coords) < 1)
                        {
                            // Get necessary information
                            const auto origin = positions(coords);
                            const auto velocity = velocities(coords);

                            // Translate to correct position
                            Eigen::Matrix<float_t, 4, 4> translation_matrix;
                            translation_matrix.setIdentity();
                            translation_matrix.col(3).head(3) = origin;

                            // Scale according to user input
                            Eigen::Matrix<float_t, 4, 4> scale_matrix;
                            scale_matrix.setIdentity();

                            switch (arrow_size)
                            {
                            case interface_deformation_glyph_aux::arrow_size_t::dynamic:
                                scale_matrix(0, 0) = scale_matrix(1, 1) = scale_matrix(2, 2) = arrow_scalar * this->timestep * velocity.norm();

                                break;
                            case interface_deformation_glyph_aux::arrow_size_t::fixed_thickness:
                                scale_matrix(0, 0) = arrow_scalar * this->timestep * velocity.norm();
                                scale_matrix(1, 1) = scale_matrix(2, 2) = arrow_fixed_scalar;

                                break;
                            case interface_deformation_glyph_aux::arrow_size_t::fixed_size:
                            default:
                                scale_matrix(0, 0) = scale_matrix(1, 1) = scale_matrix(2, 2) = arrow_fixed_scalar;
                            }

                            // Rotate to match velocity's direction
                            Eigen::Matrix<float_t, 4, 4> rotation_matrix;
                            rotation_matrix.setIdentity();

                            const auto rotation_axis = Eigen::Matrix<float_t, 3, 1>(1.0, 0.0, 0.0).cross(velocity.normalized()).normalized();
                            const auto angle = std::acos(Eigen::Matrix<float_t, 3, 1>(1.0, 0.0, 0.0).dot(velocity.normalized()));

                            auto& rotation_matrix_3x3 = rotation_matrix.block(0, 0, 3, 3);

                            Eigen::Matrix<float_t, 3, 3> cross_product_matrix;
                            cross_product_matrix <<
                                0.0, -rotation_axis[2], rotation_axis[1],
                                rotation_axis[2], 0.0, -rotation_axis[0],
                                -rotation_axis[1], rotation_axis[0], 0.0;

                            rotation_matrix_3x3 *= std::cos(angle);
                            rotation_matrix_3x3 += std::sin(angle) * cross_product_matrix;
                            rotation_matrix_3x3 += (static_cast<float_t>(1.0) - std::cos(angle)) * (rotation_axis * rotation_axis.transpose());

                            // Create object matrix
                            const math::transformer<float_t, 3> trafo(translation_matrix * rotation_matrix * scale_matrix);

                            // Instantiate glyph
                            auto glyph_instance = glyph_template->clone(trafo);

                            const auto magnitude = velocities(coords).norm();

                            #pragma omp critical
                            {
                                this->velocity_glyphs->insert(glyph_instance);

                                magnitudes->push_back(magnitude);
                            }
                        }
                    }
                }
            }

            this->velocity_glyphs->add(magnitudes, data::topology_t::OBJECT_DATA);
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_stretching_glyph_template(std::size_t circle_resolution, float_t bending,
            float_t hole_radius, const float_t offset, float_t strip_width, float_t reference_width, float_t z_offset) const -> stretching_glyph_t
        {
            // Clamp parameter values
            circle_resolution = std::max(circle_resolution, static_cast<std::size_t>(8));
            bending = std::max(bending, static_cast<float_t>(0.0));
            hole_radius = std::clamp(hole_radius, static_cast<float_t>(0.1), static_cast<float_t>(0.9));
            strip_width = std::clamp(strip_width, static_cast<float_t>(0.01), static_cast<float_t>(0.5));
            reference_width = std::clamp(reference_width, static_cast<float_t>(0.001), static_cast<float_t>(0.5));
            z_offset = std::max(z_offset, static_cast<float_t>(0.0));

            // General variables
            const std::size_t radial_resolution = 10;

            const auto circle_increment = static_cast<float_t>((2.0 * math::pi<float_t>) / circle_resolution);
            const auto radial_increment = static_cast<float_t>(1.0 - hole_radius) / radial_resolution;

            const auto disc_width = static_cast<float_t>(1.0 - hole_radius);
            const auto center_line_offset = static_cast<float_t>(0.5 * (1.0 + hole_radius));
            const auto denominator = static_cast<float_t>(std::pow(0.5 * disc_width, 2.0));

            auto bending_func = [&center_line_offset, &denominator, &bending](const float_t x) {
                return bending * (static_cast<float_t>(1.0) - std::pow(x - center_line_offset, static_cast<float_t>(2.0)) / denominator);
            };

            // Create circle-shape disc on x,y-plane
            glyph_t disc = std::make_shared<geometry::mesh<float_t>>();
            auto disc_tex_coords = std::make_shared<std::vector<Eigen::Matrix<float_t, 2, 1>>>(circle_resolution * (radial_resolution + 1));

            {
                // Create rows
                std::vector<CGAL::SM_Vertex_index> previous_indices(circle_resolution);
                std::vector<CGAL::SM_Vertex_index> current_indices(circle_resolution);

                previous_indices[0] = disc->add_point(geometry::point<float_t>(
                    hole_radius * std::cos(static_cast<float_t>(0.0)),
                    hole_radius * std::sin(static_cast<float_t>(0.0)),
                    offset + bending_func(hole_radius)));

                (*disc_tex_coords)[previous_indices[0]] <<
                    (1.0 + hole_radius * std::cos(static_cast<float_t>(0.0))) / 2.0,
                    (1.0 + hole_radius * std::sin(static_cast<float_t>(0.0))) / 2.0;

                for (std::size_t i = 1; i < circle_resolution; ++i)
                {
                    const auto angle = i * circle_increment;

                    previous_indices[i] = disc->add_point(geometry::point<float_t>(
                        hole_radius * std::cos(angle),
                        hole_radius * std::sin(angle),
                        offset + bending_func(hole_radius)));

                    (*disc_tex_coords)[previous_indices[i]] <<
                        (1.0 + hole_radius * std::cos(angle)) / 2.0,
                        (1.0 + hole_radius * std::sin(angle)) / 2.0;
                }

                for (std::size_t j = 1; j <= radial_resolution; ++j)
                {
                    const auto radius = hole_radius + j * radial_increment;

                    current_indices[0] = disc->add_point(geometry::point<float_t>(
                        radius * std::cos(static_cast<float_t>(0.0)),
                        radius * std::sin(static_cast<float_t>(0.0)),
                        offset + bending_func(radius)));

                    (*disc_tex_coords)[current_indices[0]] <<
                        (1.0 + radius * std::cos(static_cast<float_t>(0.0))) / 2.0,
                        (1.0 + radius * std::sin(static_cast<float_t>(0.0))) / 2.0;

                    for (std::size_t i = 1; i < circle_resolution; ++i)
                    {
                        const auto angle = i * circle_increment;

                        current_indices[i] = disc->add_point(geometry::point<float_t>(
                            radius * std::cos(angle),
                            radius * std::sin(angle),
                            offset + bending_func(radius)));

                        (*disc_tex_coords)[current_indices[i]] <<
                            (1.0 + radius * std::cos(angle)) / 2.0,
                            (1.0 + radius * std::sin(angle)) / 2.0;

                        disc->add_face({ previous_indices[i - 1], current_indices[i - 1], current_indices[i], previous_indices[i] });
                    }

                    disc->add_face({ previous_indices.back(), current_indices.back(), current_indices.front(), previous_indices.front() });

                    std::swap(current_indices, previous_indices);
                }
            }

            // Create strips in x- and y-direction, respectively
            glyph_t min_strip;
            glyph_t max_strip;

            auto strip_tex_coords = std::make_shared<std::vector<Eigen::Matrix<float_t, 2, 1>>>(4 * (radial_resolution + 1));
            std::fill(strip_tex_coords->begin(), strip_tex_coords->end(), Eigen::Matrix<float_t, 2, 1>(0.0, 0.0));

            {
                // Create strip in positive x-direction
                auto strip = std::make_shared<geometry::mesh<float_t>>();

                const auto y_plus = strip_width / static_cast<float_t>(2.0);
                const auto y_minus = -strip_width / static_cast<float_t>(2.0);

                // Create rows
                std::array<CGAL::SM_Vertex_index, 2> previous_indices;
                std::array<CGAL::SM_Vertex_index, 2> current_indices;

                previous_indices[0] = strip->add_point(geometry::point<float_t>(
                    hole_radius,
                    y_plus,
                    offset + z_offset + bending_func(hole_radius)));

                previous_indices[1] = strip->add_point(geometry::point<float_t>(
                    hole_radius,
                    y_minus,
                    offset + z_offset + bending_func(hole_radius)));

                for (std::size_t j = 1; j <= radial_resolution; ++j)
                {
                    const auto radius = hole_radius + j * radial_increment;

                    current_indices[0] = strip->add_point(geometry::point<float_t>(
                        radius,
                        y_plus,
                        offset + z_offset + bending_func(radius)));

                    current_indices[1] = strip->add_point(geometry::point<float_t>(
                        radius,
                        y_minus,
                        offset + z_offset + bending_func(radius)));

                    strip->add_face({ previous_indices[1], current_indices[1], current_indices[0], previous_indices[0] });

                    std::swap(current_indices, previous_indices);
                }

                // Copy strip for negative x-direction, and both y-directions
                const Eigen::Matrix<float_t, 3, 1> origin(0.0, 0.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> x_axis(1.0, 0.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> y_axis(0.0, 1.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> z_axis(0.0, 0.0, 1.0);

                math::transformer<float_t, 3> x_neg(origin, -x_axis, -y_axis, z_axis);
                math::transformer<float_t, 3> y_pos(origin, y_axis, -x_axis, z_axis);
                math::transformer<float_t, 3> y_neg(origin, -y_axis, x_axis, z_axis);

                min_strip = strip->merge(*std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(x_neg)));
                max_strip = std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(y_pos))
                    ->merge(*std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(y_neg)));
            }

            // Create a circle as reference
            glyph_t reference;

            {
                // Function to create the reference as a slim copy of the disc
                const auto ratio = reference_width / disc_width;
                const auto circle_offset = offset + z_offset + ratio * bending;

                auto slimmify = [&center_line_offset, &ratio, &circle_offset](const Eigen::Matrix<float_t, 4, 1>& vector) -> Eigen::Matrix<float_t, 4, 1>
                {
                    const auto center_line_vector = center_line_offset * vector.head<2>().normalized();
                    const auto vector_xy = center_line_vector + ratio * (vector.head<2>() - center_line_vector);

                    Eigen::Matrix<float_t, 4, 1> slim_vector;
                    slim_vector << vector_xy, circle_offset + ratio * vector[2], 1.0;

                    return slim_vector;
                };

                math::transformer<float_t, 3> trafo;
                trafo.set_preprocessing(slimmify);

                reference = std::static_pointer_cast<geometry::mesh<float_t>>(disc->clone(trafo));
            }

            return std::make_tuple(
                std::make_pair(disc, disc_tex_coords),
                std::make_pair(min_strip, strip_tex_coords),
                std::make_pair(max_strip, strip_tex_coords),
                std::make_pair(reference, disc_tex_coords));
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::instantiate_stretching_glpyh(const stretching_glyph_t& glyph_template,
            const float_t average_cell_size, float_t size_scalar, float_t scalar, float_t sharpness, const bool show_strips, const bool show_reference)
        {
            // Clamp parameter values
            size_scalar = std::max(size_scalar, static_cast<float_t>(0.00001 * average_cell_size));
            scalar = std::max(scalar, static_cast<float_t>(0.0));
            sharpness = std::max(sharpness, static_cast<float_t>(0.0));

            // Get input data
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 3>& stretching_direction_min = *this->stretching_min;
            const data::grid<float_t, float_t, 3, 3>& stretching_direction_max = *this->stretching_max;

            // Create output data arrays
            auto tex_coords_rings = std::make_shared<data::array<float_t, 2>>("Texture Coordinates");
            auto stretching_rings = std::make_shared<data::array<float_t>>("Stretching");
            auto relevance_rings = std::make_shared<data::array<float_t>>("Relevance");

            auto tex_coords_strips = std::make_shared<data::array<float_t, 2>>("Texture Coordinates");
            auto stretching_strips = std::make_shared<data::array<float_t>>("Stretching");
            auto relevance_strips = std::make_shared<data::array<float_t>>("Relevance");

            auto tex_coords_references = std::make_shared<data::array<float_t, 2>>("Texture Coordinates");
            auto stretching_references = std::make_shared<data::array<float_t>>("Stretching");
            auto relevance_references = std::make_shared<data::array<float_t>>("Relevance");

            // Iterate over all interface cells
            #pragma omp parallel for schedule(dynamic) default(none) shared(size_scalar, scalar)
            for (long long z = static_cast<long long>(vof.get_extent()[2].first);
                z <= static_cast<long long>(vof.get_extent()[2].second); ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (vof(coords) > 0 && vof(coords) < 1)
                        {
                            // Get necessary information
                            const auto origin = positions(coords);
                            const auto normal = -gradients(coords).normalized();
                            const auto stretching_min = stretching_direction_min(coords);
                            const auto stretching_max = stretching_direction_max(coords);

                            const auto min = stretching_min.norm();
                            const auto max = stretching_max.norm();

                            // Rotate into basis defined by the interface normal and the eigenvectors
                            const math::transformer<float_t, 3> translation_and_rotation(origin,
                                stretching_min.normalized(), stretching_max.normalized(), normal);

                            // Calculate and apply anisotropy for superellipse
                            const auto anisotropy = (max - min) / (max + min);
                            const auto anisotropy_exp = std::pow(1.0 - anisotropy, sharpness);

                            auto superellipse = [anisotropy_exp](const Eigen::Matrix<float_t, 4, 1>& vector) -> Eigen::Matrix<float_t, 4, 1>
                            {
                                const auto radius = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);

                                Eigen::Matrix<float_t, 4, 1> superellipse_vector;

                                superellipse_vector <<
                                    radius * std::copysign(std::pow(std::abs(vector[0]) / radius, anisotropy_exp), vector[0]),
                                    radius * std::copysign(std::pow(std::abs(vector[1]) / radius, anisotropy_exp), vector[1]),
                                    vector[2],
                                    vector[3];

                                return superellipse_vector;
                            };

                            // Scale
                            Eigen::Matrix<float_t, 4, 4> scale_matrix, reference_scale_matrix;
                            scale_matrix.setIdentity();
                            reference_scale_matrix.setIdentity();

                            scale_matrix(0, 0) = size_scalar * average_cell_size * min * ((min > 1.0) ? scalar : (1.0 / scalar));
                            scale_matrix(1, 1) = size_scalar * average_cell_size * max * ((max > 1.0) ? scalar : (1.0 / scalar));
                            scale_matrix(2, 2) = size_scalar * average_cell_size;

                            reference_scale_matrix(0, 0) = size_scalar * average_cell_size;
                            reference_scale_matrix(1, 1) = size_scalar * average_cell_size;
                            reference_scale_matrix(2, 2) = size_scalar * average_cell_size;

                            // Create object matrix
                            math::transformer<float_t, 3> trafo(translation_and_rotation * math::transformer<float_t, 3>(scale_matrix));
                            const math::transformer<float_t, 3> trafo_ref(translation_and_rotation * math::transformer<float_t, 3>(reference_scale_matrix));

                            trafo.set_preprocessing(superellipse);

                            // Instantiate glyph
                            auto stretching_glyph_ring_instance = std::get<0>(glyph_template).first->clone(trafo);

                            std::shared_ptr<geometry::geometric_object<float_t>>
                                stretching_glyph_strip_1_instance, stretching_glyph_strip_2_instance,
                                stretching_glyph_reference_1_instance, stretching_glyph_reference_2_instance;

                            if (show_strips)
                            {
                                stretching_glyph_strip_1_instance = std::get<1>(glyph_template).first->clone(trafo);
                                stretching_glyph_strip_2_instance = std::get<2>(glyph_template).first->clone(trafo);
                            }

                            if (show_reference)
                            {
                                stretching_glyph_reference_1_instance = std::get<3>(glyph_template).first->clone(trafo_ref);
                                stretching_glyph_reference_2_instance = std::get<3>(glyph_template).first->clone(trafo);
                            }

                            const auto min_abs_exp = (min < 1.0) ? (static_cast<float_t>(1.0) / min) : min;
                            const auto max_abs_exp = (max < 1.0) ? (static_cast<float_t>(1.0) / max) : max;

                            #pragma omp critical
                            {
                                this->stretching_glyph_rings->insert(stretching_glyph_ring_instance);

                                tex_coords_rings->insert_back(std::get<0>(glyph_template).second->begin(), std::get<0>(glyph_template).second->end());
                                stretching_rings->push_back(1.0);
                                relevance_rings->push_back(min_abs_exp * max_abs_exp);

                                if (show_strips)
                                {
                                    this->stretching_glyph_strips->insert(stretching_glyph_strip_1_instance);
                                    this->stretching_glyph_strips->insert(stretching_glyph_strip_2_instance);

                                    tex_coords_strips->insert_back(std::get<1>(glyph_template).second->begin(), std::get<1>(glyph_template).second->end());
                                    tex_coords_strips->insert_back(std::get<2>(glyph_template).second->begin(), std::get<2>(glyph_template).second->end());

                                    stretching_strips->push_back(min);
                                    stretching_strips->push_back(max);

                                    relevance_strips->push_back(min_abs_exp * max_abs_exp);
                                    relevance_strips->push_back(min_abs_exp * max_abs_exp);
                                }

                                if (show_reference)
                                {
                                    this->stretching_glyph_references->insert(stretching_glyph_reference_1_instance);
                                    this->stretching_glyph_references->insert(stretching_glyph_reference_2_instance);

                                    tex_coords_references->insert_back(std::get<3>(glyph_template).second->begin(), std::get<3>(glyph_template).second->end());
                                    tex_coords_references->insert_back(std::get<3>(glyph_template).second->begin(), std::get<3>(glyph_template).second->end());

                                    stretching_references->push_back(std::numeric_limits<float_t>::infinity());     // TODO
                                    stretching_references->push_back(std::numeric_limits<float_t>::quiet_NaN());    // TODO

                                    relevance_references->push_back(min_abs_exp * max_abs_exp);
                                    relevance_references->push_back(min_abs_exp * max_abs_exp);
                                }
                            }
                        }
                    }
                }
            }

            this->stretching_glyph_rings->add(tex_coords_rings, data::topology_t::TEXTURE_COORDINATES);
            this->stretching_glyph_rings->add(stretching_rings, data::topology_t::OBJECT_DATA);
            this->stretching_glyph_rings->add(relevance_rings, data::topology_t::OBJECT_DATA);

            this->stretching_glyph_strips->add(tex_coords_strips, data::topology_t::TEXTURE_COORDINATES);
            this->stretching_glyph_strips->add(stretching_strips, data::topology_t::OBJECT_DATA);
            this->stretching_glyph_strips->add(relevance_strips, data::topology_t::OBJECT_DATA);

            this->stretching_glyph_references->add(tex_coords_references, data::topology_t::TEXTURE_COORDINATES);
            this->stretching_glyph_references->add(stretching_references, data::topology_t::OBJECT_DATA);
            this->stretching_glyph_references->add(relevance_references, data::topology_t::OBJECT_DATA);
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_bending_glyph_template(std::size_t circle_resolution,
            std::size_t polynomial_resolution, const float_t offset, float_t strip_width, float_t z_offset) const -> bending_glyph_t
        {
            // Clamp parameter values
            circle_resolution = std::max(circle_resolution, static_cast<std::size_t>(8));
            polynomial_resolution = std::max(polynomial_resolution, static_cast<std::size_t>(1));
            strip_width = std::clamp(strip_width, static_cast<float_t>(0.01), static_cast<float_t>(0.5));
            z_offset = std::max(z_offset, static_cast<float_t>(0.0));

            // Create circle-shape disc on x,y-plane
            glyph_t disc = std::make_shared<geometry::mesh<float_t>>();

            {
                const auto circle_increment = static_cast<float_t>((2.0 * math::pi<float_t>) / circle_resolution);
                const auto polynomial_increment = static_cast<float_t>(1.0 / polynomial_resolution);

                // Create inner row
                std::vector<CGAL::SM_Vertex_index> previous_indices(circle_resolution);

                {
                    const auto center = disc->add_point(geometry::point<float_t>(0.0, 0.0, offset + z_offset));

                    previous_indices[0] = disc->add_point(geometry::point<float_t>(
                        polynomial_increment * std::cos(static_cast<float_t>(0.0)),
                        polynomial_increment * std::sin(static_cast<float_t>(0.0)),
                        offset + z_offset));

                    for (std::size_t i = 1; i < circle_resolution; ++i)
                    {
                        const auto angle = i * circle_increment;

                        previous_indices[i] = disc->add_point(geometry::point<float_t>(
                            polynomial_increment * std::cos(angle),
                            polynomial_increment * std::sin(angle),
                            offset + z_offset));

                        disc->add_face(center, previous_indices[i - 1], previous_indices[i]);
                    }

                    disc->add_face(center, previous_indices.back(), previous_indices.front());
                }

                // Create outer rows
                if (polynomial_resolution > 1)
                {
                    std::vector<CGAL::SM_Vertex_index> current_indices(circle_resolution);

                    for (std::size_t j = 2; j <= polynomial_resolution; ++j)
                    {
                        auto radius = j * polynomial_increment;

                        current_indices[0] = disc->add_point(geometry::point<float_t>(
                            radius * std::cos(static_cast<float_t>(0.0)),
                            radius * std::sin(static_cast<float_t>(0.0)),
                            offset + z_offset));

                        for (std::size_t i = 1; i < circle_resolution; ++i)
                        {
                            const auto angle = i * circle_increment;

                            current_indices[i] = disc->add_point(geometry::point<float_t>(
                                radius * std::cos(angle),
                                radius * std::sin(angle),
                                offset + z_offset));

                            disc->add_face({ previous_indices[i - 1], current_indices[i - 1], current_indices[i], previous_indices[i] });
                        }

                        disc->add_face({ previous_indices.back(), current_indices.back(), current_indices.front(), previous_indices.front() });

                        std::swap(current_indices, previous_indices);
                    }
                }
            }

            // Create strips in x- and y-direction, respectively
            glyph_t min_strip;
            glyph_t max_strip;

            {
                // Create strip in positive x-direction
                auto strip = std::make_shared<geometry::mesh<float_t>>();

                const auto increment = static_cast<float_t>((1.0 - strip_width) / polynomial_resolution);

                const auto y_plus = strip_width / static_cast<float_t>(2.0);
                const auto y_minus = -strip_width / static_cast<float_t>(2.0);

                auto index_prev_1 = strip->add_point(geometry::point<float_t>(strip_width, y_plus, offset + 2.0 * z_offset));
                auto index_prev_2 = strip->add_point(geometry::point<float_t>(strip_width, y_minus, offset + 2.0 * z_offset));

                for (std::size_t i = 1; i <= polynomial_resolution; ++i)
                {
                    const auto x = i * increment + strip_width;

                    const auto index_1 = strip->add_point(geometry::point<float_t>(x, y_plus, offset + 2.0 * z_offset));
                    const auto index_2 = strip->add_point(geometry::point<float_t>(x, y_minus, offset + 2.0 * z_offset));

                    strip->add_face({ index_prev_1, index_prev_2, index_2, index_1 });

                    index_prev_1 = index_1;
                    index_prev_2 = index_2;
                }

                // Copy strip for negative x-direction, and both y-directions
                const Eigen::Matrix<float_t, 3, 1> origin(0.0, 0.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> x_axis(1.0, 0.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> y_axis(0.0, 1.0, 0.0);
                const Eigen::Matrix<float_t, 3, 1> z_axis(0.0, 0.0, 1.0);

                math::transformer<float_t, 3> x_neg(origin, -x_axis, -y_axis, z_axis);
                math::transformer<float_t, 3> y_pos(origin, y_axis, -x_axis, z_axis);
                math::transformer<float_t, 3> y_neg(origin, -y_axis, x_axis, z_axis);

                min_strip = strip->merge(*std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(x_neg)));
                max_strip = std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(y_pos))
                    ->merge(*std::static_pointer_cast<geometry::mesh<float_t>>(strip->clone(y_neg)));
            }

            // Return glyphs
            return std::make_tuple(disc, min_strip, max_strip);
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::instantiate_bending_glyph(const bending_glyph_t& glyph_template,
            const float_t average_cell_size, float_t size_scalar, float_t scalar, const bool show_strips)
        {
            // Clamp parameter values
            size_scalar = std::max(size_scalar, static_cast<float_t>(0.00001 * average_cell_size));
            scalar = std::max(scalar, static_cast<float_t>(0.00001 * average_cell_size));

            // Get input data
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 1>& bending_min = *this->bending_min;
            const data::grid<float_t, float_t, 3, 1>& bending_max = *this->bending_max;
            const data::grid<float_t, float_t, 3, 3>& bending_direction_min = *this->bending_direction_min;
            const data::grid<float_t, float_t, 3, 3>& bending_direction_max = *this->bending_direction_max;
            const data::grid<float_t, float_t, 3, 3>& bending_polynomial = *this->bending_polynomial;

            // Create output data arrays
            auto bending_discs = std::make_shared<data::array<float_t>>("Bending");
            auto relevance_discs = std::make_shared<data::array<float_t>>("Relevance");

            auto bending_strips = std::make_shared<data::array<float_t>>("Bending");
            auto relevance_strips = std::make_shared<data::array<float_t>>("Relevance");

            // Iterate over all interface cells
            #pragma omp parallel for schedule(dynamic) default(none) shared(size_scalar, scalar)
            for (long long z = static_cast<long long>(vof.get_extent()[2].first);
                z <= static_cast<long long>(vof.get_extent()[2].second); ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (vof(coords) > 0 && vof(coords) < 1)
                        {
                            // Get necessary information
                            const auto origin = positions(coords);
                            const auto normal = -gradients(coords).normalized();
                            const auto min = bending_min(coords);
                            const auto max = bending_max(coords);
                            const auto direction_min = bending_direction_min(coords).normalized();
                            const auto direction_max = bending_direction_max(coords).normalized();
                            const auto polynomial = bending_polynomial(coords);

                            // Rotate into basis defined by the interface normal and the eigenvectors
                            const math::transformer<float_t, 3> translation_and_rotation(origin, direction_min, direction_max, normal);

                            // Scale
                            Eigen::Matrix<float_t, 4, 4> scale_matrix;
                            scale_matrix.setIdentity();

                            scale_matrix(0, 0) = size_scalar * average_cell_size;
                            scale_matrix(1, 1) = size_scalar * average_cell_size;
                            scale_matrix(2, 2) = size_scalar * average_cell_size;

                            // Deform, such that the glyph indicates bending
                            auto bend = [&polynomial, scalar](const Eigen::Matrix<float_t, 4, 1>& vector) -> Eigen::Matrix<float_t, 4, 1>
                            {
                                // If the mixed part of the polynomial (a_xy) is not zero, the eigenvectors do not align with the
                                // x- and y-axis. This "rotation" has to be removed before sampling the polynomial.
                                math::transformer<float_t, 4> to_origin = math::transformer<float_t, 4>::unit();

                                if (!math::equals(polynomial[0], static_cast<float_t>(0.0L)))
                                {
                                    const auto meanCurv = polynomial[1] + polynomial[2];
                                    const auto gaussCurv = static_cast<float_t>(4.0L) * polynomial[1] * polynomial[2] - polynomial[0] * polynomial[0];

                                    const auto first_curvature = meanCurv - std::sqrt(meanCurv * meanCurv - gaussCurv);
                                    const auto second_curvature = meanCurv + std::sqrt(meanCurv * meanCurv - gaussCurv);

                                    Eigen::Matrix<float_t, 3, 1> x_axis, y_axis;
                                    x_axis << polynomial[0], std::min(first_curvature, second_curvature) - static_cast<float_t>(2.0L) * polynomial[1], 0.0;
                                    y_axis << polynomial[0], std::max(first_curvature, second_curvature) - static_cast<float_t>(2.0L) * polynomial[1], 0.0;

                                    to_origin = math::transformer<float_t, 4>(Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0),
                                        x_axis.normalized(), y_axis.normalized(), Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0));
                                }

                                const Eigen::Matrix<float_t, 4, 1> rotated_vector = to_origin.transform(vector);

                                Eigen::Matrix<float_t, 4, 1> bent_vector;
                                bent_vector << rotated_vector[0], rotated_vector[1],
                                    rotated_vector[2] + scalar * (polynomial[0] * rotated_vector[0] * rotated_vector[1] + polynomial[1] * rotated_vector[0] * rotated_vector[0]
                                        + polynomial[2] * rotated_vector[1] * rotated_vector[1]) + rotated_vector[2],
                                    rotated_vector[3];

                                return to_origin.transform_inverse(bent_vector);
                            };

                            // Create object matrix
                            math::transformer<float_t, 3> trafo(translation_and_rotation * math::transformer<float_t, 3>(scale_matrix));
                            trafo.set_preprocessing(bend);

                            // Instantiate glyph
                            auto bending_glyph_disc_instance = std::get<0>(glyph_template)->clone(trafo);

                            std::shared_ptr<geometry::geometric_object<float_t>>
                                bending_glyph_strip_1_instance, bending_glyph_strip_2_instance;

                            if (show_strips)
                            {
                                bending_glyph_strip_1_instance = std::get<1>(glyph_template)->clone(trafo);
                                bending_glyph_strip_2_instance = std::get<2>(glyph_template)->clone(trafo);
                            }

                            #pragma omp critical
                            {
                                this->bending_glyph_discs->insert(bending_glyph_disc_instance);

                                bending_discs->push_back(0.0);
                                relevance_discs->push_back(std::abs(min) + std::abs(max));

                                if (show_strips)
                                {
                                    this->bending_glyph_strips->insert(bending_glyph_strip_1_instance);
                                    this->bending_glyph_strips->insert(bending_glyph_strip_2_instance);

                                    bending_strips->push_back(min);
                                    bending_strips->push_back(max);

                                    relevance_strips->push_back(std::abs(min) + std::abs(max));
                                    relevance_strips->push_back(std::abs(min) + std::abs(max));
                                }
                            }
                        }
                    }
                }
            }

            this->bending_glyph_discs->add(bending_discs, data::topology_t::OBJECT_DATA);
            this->bending_glyph_discs->add(relevance_discs, data::topology_t::OBJECT_DATA);

            this->bending_glyph_strips->add(bending_strips, data::topology_t::OBJECT_DATA);
            this->bending_glyph_strips->add(relevance_strips, data::topology_t::OBJECT_DATA);
        }
    }
}
