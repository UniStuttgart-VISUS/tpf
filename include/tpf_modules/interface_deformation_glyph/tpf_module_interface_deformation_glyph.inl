#include "tpf_module_interface_deformation_glyph.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_mesh.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_math_defines.h"
#include "tpf/math/tpf_transformer.h"
#include "tpf/math/tpf_vector.h"

#include <cmath>
#include <memory>
#include <string>
#include <tuple>
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
            data::polydata<float_t>& stretching_glyphs, data::polydata<float_t>& bending_glyphs)
        {
            this->velocity_glyphs = &velocity_glyphs;
            this->stretching_glyphs = &stretching_glyphs;
            this->bending_glyphs = &bending_glyphs;
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::set_algorithm_parameters(const bool velocity_glyph,
            const bool stretching_glyph, const bool bending_glyph, const float_t timestep,
            const interface_deformation_glyph_aux::arrow_size_t arrow_size, const float_t arrow_scalar,
            const float_t arrow_fixed_scalar, const int arrow_resolution, const float_t arrow_ratio, const float_t arrow_thickness,
            const int bending_disc_resolution, const int bending_polygonal_resolution, const float_t bending_strip_size,
            const float_t bending_size_scalar, const float_t bending_scalar)
        {
            this->velocity_glyph = velocity_glyph;
            this->stretching_glyph = stretching_glyph;
            this->bending_glyph = bending_glyph;

            this->timestep = timestep;

            this->arrow_size = arrow_size;
            this->arrow_scalar = arrow_scalar;
            this->arrow_fixed_scalar = arrow_fixed_scalar;
            this->arrow_resolution = arrow_resolution;
            this->arrow_ratio = arrow_ratio;
            this->arrow_thickness = arrow_thickness;

            this->bending_disc_resolution = bending_disc_resolution;
            this->bending_polygonal_resolution = bending_polygonal_resolution;
            this->bending_strip_size = bending_strip_size;
            this->bending_size_scalar = bending_size_scalar;
            this->bending_scalar = bending_scalar;
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
                instantiate_velocity_glyphs(create_velocity_glyph_template(this->arrow_resolution, this->arrow_ratio, this->arrow_thickness),
                    this->arrow_size, this->arrow_scalar, this->arrow_fixed_scalar);
            }

            // Create stretching glyph
            if (this->stretching_glyph && this->gradients != nullptr && this->stretching_min != nullptr && this->stretching_max != nullptr)
            {
                // TODO

                // Stretch according to eigenvalues in x,y-direction
                // Rotate such that the new basis are the eigenvectors and the interface normal
                //   parameter: average_cell_size
                //              stretch_min, stretch_max (relative to average_cell_size)
            }

            // Create bending glyph
            if (this->bending_glyph && this->gradients != nullptr && this->bending_min != nullptr && this->bending_max != nullptr &&
                this->bending_direction_min != nullptr && this->bending_direction_max != nullptr)
            {
                instantiate_bending_glyph(create_bending_glyph_template(this->bending_disc_resolution,
                    this->bending_polygonal_resolution, this->bending_strip_size),
                    average_cell_size, this->bending_size_scalar, this->bending_scalar);
            }
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_velocity_glyph_template(const std::size_t resolution,
            const float_t shaft_tip_ratio, const float_t thickness_ratio) const -> velocity_glyph_t
        {
            const auto shaft_thickness = thickness_ratio;
            const auto tip_thickness = static_cast<float_t>(2.0) * thickness_ratio;

            const auto increment = static_cast<float_t>((2.0 * math::pi<float_t>) / resolution);

            // Create arrow glyph...
            velocity_glyph_t arrow_glyph = std::make_shared<geometry::mesh<float_t>>();

            {
                // ... shaft
                {
                    std::vector<CGAL::SM_Vertex_index> origin_indices(resolution);
                    std::vector<CGAL::SM_Vertex_index> target_indices(resolution);

                    const auto cap_center = arrow_glyph->add_point(geometry::point<float_t>(0.0, 0.0, 0.0));

                    origin_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        0.0,
                        shaft_thickness * std::cos(static_cast<float_t>(0.0)),
                        shaft_thickness * std::sin(static_cast<float_t>(0.0))));

                    target_indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        shaft_tip_ratio,
                        shaft_thickness * std::cos(static_cast<float_t>(0.0)),
                        shaft_thickness * std::sin(static_cast<float_t>(0.0))));

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        origin_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            0.0,
                            shaft_thickness * std::cos(angle),
                            shaft_thickness * std::sin(angle)));

                        target_indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            shaft_tip_ratio,
                            shaft_thickness * std::cos(angle),
                            shaft_thickness * std::sin(angle)));

                        arrow_glyph->add_face(cap_center, origin_indices[i - 1], origin_indices[i]);
                        arrow_glyph->add_face({ origin_indices[i - 1], target_indices[i - 1], target_indices[i], origin_indices[i] });
                    }

                    arrow_glyph->add_face(cap_center, origin_indices.back(), origin_indices.front());
                    arrow_glyph->add_face({ origin_indices.back(), target_indices.back(), target_indices.front(), origin_indices.front() });
                }

                // ... tip
                {
                    std::vector<CGAL::SM_Vertex_index> indices(resolution);

                    const auto tip_center = arrow_glyph->add_point(geometry::point<float_t>(1.0, 0.0, 0.0));
                    const auto tip_cap_center = arrow_glyph->add_point(geometry::point<float_t>(shaft_tip_ratio, 0.0, 0.0));

                    indices[0] = arrow_glyph->add_point(geometry::point<float_t>(
                        shaft_tip_ratio,
                        tip_thickness * std::cos(static_cast<float_t>(0.0)),
                        tip_thickness * std::sin(static_cast<float_t>(0.0))));

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        indices[i] = arrow_glyph->add_point(geometry::point<float_t>(
                            shaft_tip_ratio,
                            tip_thickness * std::cos(angle),
                            tip_thickness * std::sin(angle)));

                        arrow_glyph->add_face(tip_cap_center, indices[i - 1], indices[i]);
                        arrow_glyph->add_face(tip_center, indices[i], indices[i - 1]);
                    }

                    arrow_glyph->add_face(tip_cap_center, indices.back(), indices.front());
                    arrow_glyph->add_face(tip_center, indices.front(), indices.back());
                }
            }

            return arrow_glyph;
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::instantiate_velocity_glyphs(const velocity_glyph_t& glyph_template,
            const interface_deformation_glyph_aux::arrow_size_t arrow_size, const float_t arrow_scalar,
            const float_t arrow_fixed_scalar)
        {
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& velocities = *this->velocities;

            // Iterate over all interface cells
            std::size_t num_interface_cells = 0;

            for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
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
                            this->velocity_glyphs->insert(glyph_template->clone(trafo));

                            ++num_interface_cells;
                        }
                    }
                }
            }

            // Create data array
            auto magnitudes = std::make_shared<data::array<float_t>>("Velocity Magnitude");
            magnitudes->reserve(num_interface_cells);

            if (num_interface_cells != 0)
            {
                for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
                {
                    for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                    {
                        for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            if (vof(coords) > 0 && vof(coords) < 1)
                            {
                                const auto magnitude = velocities(coords).norm();

                                magnitudes->push_back(magnitude);
                            }
                        }
                    }
                }
            }

            this->velocity_glyphs->add(magnitudes, data::topology_t::OBJECT_DATA);
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_bending_glyph_template(const std::size_t circle_resolution,
            const std::size_t polygonal_resolution, const float_t strip_size) const -> bending_glyph_t
        {
            // Create circle-shape disc on x,y-plane
            glyph_t disc = std::make_shared<geometry::mesh<float_t>>();

            {
                const auto circle_increment = static_cast<float_t>((2.0 * math::pi<float_t>) / circle_resolution);
                const auto polygonal_increment = static_cast<float_t>(1.0 / polygonal_resolution);

                // Create inner row
                std::vector<CGAL::SM_Vertex_index> previous_indices(circle_resolution);

                {
                    const auto center = disc->add_point(geometry::point<float_t>(0.0, 0.0, 0.0));

                    previous_indices[0] = disc->add_point(geometry::point<float_t>(
                        polygonal_increment * std::cos(static_cast<float_t>(0.0)),
                        polygonal_increment * std::sin(static_cast<float_t>(0.0)),
                        0.0));

                    for (std::size_t i = 1; i < circle_resolution; ++i)
                    {
                        const auto angle = i * circle_increment;

                        previous_indices[i] = disc->add_point(geometry::point<float_t>(
                            polygonal_increment * std::cos(angle),
                            polygonal_increment * std::sin(angle),
                            0.0));

                        disc->add_face(center, previous_indices[i - 1], previous_indices[i]);
                    }

                    disc->add_face(center, previous_indices.back(), previous_indices.front());
                }

                // Create outer rows
                if (polygonal_resolution > 1)
                {
                    std::vector<CGAL::SM_Vertex_index> current_indices(circle_resolution);

                    for (std::size_t j = 2; j <= polygonal_resolution; ++j)
                    {
                        auto radius = j * polygonal_increment;

                        current_indices[0] = disc->add_point(geometry::point<float_t>(
                            radius * std::cos(static_cast<float_t>(0.0)),
                            radius * std::sin(static_cast<float_t>(0.0)),
                            0.0));

                        for (std::size_t i = 1; i < circle_resolution; ++i)
                        {
                            const auto angle = i * circle_increment;

                            current_indices[i] = disc->add_point(geometry::point<float_t>(
                                radius * std::cos(angle),
                                radius * std::sin(angle),
                                0.0));

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

                const auto increment = static_cast<float_t>((1.0 - strip_size) / polygonal_resolution);

                const auto y_plus = strip_size / static_cast<float_t>(2.0);
                const auto y_minus = -strip_size / static_cast<float_t>(2.0);
                const auto z_offset = static_cast<float_t>(0.02);

                auto index_prev_1 = strip->add_point(geometry::point<float_t>(strip_size, y_plus, z_offset));
                auto index_prev_2 = strip->add_point(geometry::point<float_t>(strip_size, y_minus, z_offset));

                for (std::size_t i = 1; i <= polygonal_resolution; ++i)
                {
                    const auto x = i * increment + strip_size;

                    const auto index_1 = strip->add_point(geometry::point<float_t>(x, y_plus, z_offset));
                    const auto index_2 = strip->add_point(geometry::point<float_t>(x, y_minus, z_offset));

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
            const float_t average_cell_size, const float_t size_scalar, const float_t scalar)
        {
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 1>& bending_min = *this->bending_min;
            const data::grid<float_t, float_t, 3, 1>& bending_max = *this->bending_max;
            const data::grid<float_t, float_t, 3, 3>& bending_direction_min = *this->bending_direction_min;
            const data::grid<float_t, float_t, 3, 3>& bending_direction_max = *this->bending_direction_max;
            const data::grid<float_t, float_t, 3, 3>& bending_polynomial = *this->bending_polynomial;

            // Iterate over all interface cells
            std::size_t num_interface_cells = 0;

            for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
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
                                    scalar * (polynomial[0] * rotated_vector[0] * rotated_vector[1] + polynomial[1] * rotated_vector[0] * rotated_vector[0]
                                        + polynomial[2] * rotated_vector[1] * rotated_vector[1]) + rotated_vector[2],
                                    rotated_vector[3];

                                return to_origin.transform_inverse(bent_vector);
                            };

                            // Create object matrix
                            math::transformer<float_t, 3> trafo(translation_and_rotation * math::transformer<float_t, 3>(scale_matrix));
                            trafo.set_preprocessing(bend);

                            // Instantiate glyph
                            this->bending_glyphs->insert(std::vector<std::shared_ptr<geometry::geometric_object<float_t>>>{
                                std::get<0>(glyph_template)->clone(trafo),
                                std::get<1>(glyph_template)->clone(trafo),
                                std::get<2>(glyph_template)->clone(trafo) });

                            ++num_interface_cells;
                        }
                    }
                }
            }

            // Create data array
            auto bending = std::make_shared<data::array<float_t>>("Bending");
            bending->reserve(num_interface_cells * 3);

            if (num_interface_cells != 0)
            {
                for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second; ++z)
                {
                    for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second; ++y)
                    {
                        for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            if (vof(coords) > 0 && vof(coords) < 1)
                            {
                                const auto min = bending_min(coords);
                                const auto max = bending_max(coords);

                                bending->push_back((min + max) / 2.0);
                                bending->push_back(min);
                                bending->push_back(max);
                            }
                        }
                    }
                }
            }

            this->bending_glyphs->add(bending, data::topology_t::OBJECT_DATA);
        }
    }
}
