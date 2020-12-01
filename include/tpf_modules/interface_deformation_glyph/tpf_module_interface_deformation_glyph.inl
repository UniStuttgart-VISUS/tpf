#include "tpf_module_interface_deformation_glyph.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_math_defines.h"
#include "tpf/math/tpf_transformer.h"
#include "tpf/math/tpf_vector.h"

#include <algorithm>
#include <cmath>
#include <string>

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
            velocity_glyph_t arrow_glyph;

            {
                // ... shaft
                std::vector<std::shared_ptr<geometry::triangle<float_t>>> shaft(resolution * 2);
                std::vector<std::shared_ptr<geometry::triangle<float_t>>> shaft_cap(resolution);

                {
                    auto y_prev = shaft_thickness * std::cos(static_cast<float_t>(0.0));
                    auto z_prev = shaft_thickness * std::sin(static_cast<float_t>(0.0));

                    const geometry::point<float_t> cap_center(0.0, 0.0, 0.0);

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        const auto y = shaft_thickness * std::cos(angle);
                        const auto z = shaft_thickness * std::sin(angle);

                        const geometry::point<float_t> p1(0.0, y_prev, z_prev);
                        const geometry::point<float_t> p2(shaft_tip_ratio, y_prev, z_prev);
                        const geometry::point<float_t> p3(0.0, y, z);
                        const geometry::point<float_t> p4(shaft_tip_ratio, y, z);

                        shaft[i * 2 + 0] = std::make_shared<geometry::triangle<float_t>>(p1, p3, p2);
                        shaft[i * 2 + 1] = std::make_shared<geometry::triangle<float_t>>(p3, p4, p2);

                        shaft_cap[i] = std::make_shared<geometry::triangle<float_t>>(p1, p3, cap_center);

                        y_prev = y;
                        z_prev = z;
                    }

                    {
                        const auto y = shaft_thickness * std::cos(static_cast<float_t>(0.0));
                        const auto z = shaft_thickness * std::sin(static_cast<float_t>(0.0));

                        const geometry::point<float_t> p1(0.0, y_prev, z_prev);
                        const geometry::point<float_t> p2(shaft_tip_ratio, y_prev, z_prev);
                        const geometry::point<float_t> p3(0.0, y, z);
                        const geometry::point<float_t> p4(shaft_tip_ratio, y, z);

                        shaft[0] = std::make_shared<geometry::triangle<float_t>>(p1, p3, p2);
                        shaft[1] = std::make_shared<geometry::triangle<float_t>>(p3, p4, p2);

                        shaft_cap[0] = std::make_shared<geometry::triangle<float_t>>(p1, p3, cap_center);
                    }
                }

                // ... tip
                std::vector<std::shared_ptr<geometry::triangle<float_t>>> tip(resolution);
                std::vector<std::shared_ptr<geometry::triangle<float_t>>> tip_cap(resolution);

                {
                    auto y_prev = tip_thickness * std::cos(static_cast<float_t>(0.0));
                    auto z_prev = tip_thickness * std::sin(static_cast<float_t>(0.0));

                    const geometry::point<float_t> tip_point(1.0, 0.0, 0.0);
                    const geometry::point<float_t> cap_center(shaft_tip_ratio, 0.0, 0.0);

                    for (std::size_t i = 1; i < resolution; ++i)
                    {
                        const auto angle = i * increment;

                        const auto y = tip_thickness * std::cos(angle);
                        const auto z = tip_thickness * std::sin(angle);

                        const geometry::point<float_t> p1(shaft_tip_ratio, y_prev, z_prev);
                        const geometry::point<float_t> p2(shaft_tip_ratio, y, z);

                        tip[i] = std::make_shared<geometry::triangle<float_t>>(p1, p2, tip_point);
                        tip_cap[i] = std::make_shared<geometry::triangle<float_t>>(p1, p2, cap_center);

                        y_prev = y;
                        z_prev = z;
                    }

                    {
                        const auto y = tip_thickness * std::cos(static_cast<float_t>(0.0));
                        const auto z = tip_thickness * std::sin(static_cast<float_t>(0.0));

                        const geometry::point<float_t> p1(shaft_tip_ratio, y_prev, z_prev);
                        const geometry::point<float_t> p2(shaft_tip_ratio, y, z);

                        tip[0] = std::make_shared<geometry::triangle<float_t>>(p1, p2, tip_point);
                        tip_cap[0] = std::make_shared<geometry::triangle<float_t>>(p1, p2, cap_center);
                    }
                }

                // Create single vector
                arrow_glyph.reserve(shaft.size() + tip.size() + 2);

                arrow_glyph.insert(arrow_glyph.end(), shaft.begin(), shaft.end());
                arrow_glyph.insert(arrow_glyph.end(), shaft_cap.begin(), shaft_cap.end());
                arrow_glyph.insert(arrow_glyph.end(), tip.begin(), tip.end());
                arrow_glyph.insert(arrow_glyph.end(), tip_cap.begin(), tip_cap.end());
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

            velocity_glyph_t instance(glyph_template.size());

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
                            std::transform(glyph_template.begin(), glyph_template.end(), instance.begin(),
                                [&trafo](const std::shared_ptr<geometry::geometric_object<float_t>> obj) { return obj->clone(trafo); });

                            this->velocity_glyphs->insert(instance);

                            ++num_interface_cells;
                        }
                    }
                }
            }

            // Create data array
            auto magnitudes = std::make_shared<tpf::data::array<float_t>>("Velocity Magnitude");
            magnitudes->reserve(num_interface_cells * glyph_template.size());

            for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second && num_interface_cells != 0; ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second && num_interface_cells != 0; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second && num_interface_cells != 0; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (vof(coords) > 0 && vof(coords) < 1)
                        {
                            const auto magnitude = velocities(coords).norm();

                            for (std::size_t i = 0; i < glyph_template.size(); ++i)
                            {
                                magnitudes->push_back(magnitude);
                            }
                        }
                    }
                }
            }

            this->velocity_glyphs->add(magnitudes, tpf::data::topology_t::OBJECT_DATA);
        }

        template <typename float_t>
        inline auto interface_deformation_glyph<float_t>::create_bending_glyph_template(const std::size_t circle_resolution,
            const std::size_t polygonal_resolution, const float_t strip_size) const -> bending_glyph_t
        {
            using single_glyph_t = std::vector<std::shared_ptr<geometry::triangle<float_t>>>;

            // Create circle-shape disc on x,y-plane
            glyph_t disc;
            disc.reserve(circle_resolution + 2 * circle_resolution * (polygonal_resolution - 1));

            {
                // Create inner row
                single_glyph_t inner(circle_resolution);

                {
                    const auto circle_increment = static_cast<float_t>((2.0 * math::pi<float_t>) / circle_resolution);
                    const auto polygonal_increment = static_cast<float_t>(1.0 / polygonal_resolution);

                    auto x_prev = polygonal_increment * std::cos(static_cast<float_t>(0.0));
                    auto y_prev = polygonal_increment * std::sin(static_cast<float_t>(0.0));

                    const geometry::point<float_t> center(0.0, 0.0, 0.0);

                    for (std::size_t i = 1; i < circle_resolution; ++i)
                    {
                        const auto angle = i * circle_increment;

                        const auto x = polygonal_increment * std::cos(angle);
                        const auto y = polygonal_increment * std::sin(angle);

                        const geometry::point<float_t> p1(x_prev, y_prev, 0.0);
                        const geometry::point<float_t> p2(x, y, 0.0);

                        inner[i] = std::make_shared<geometry::triangle<float_t>>(center, p1, p2);

                        x_prev = x;
                        y_prev = y;
                    }

                    const auto x = polygonal_increment * std::cos(static_cast<float_t>(0.0));
                    const auto y = polygonal_increment * std::sin(static_cast<float_t>(0.0));

                    const geometry::point<float_t> p1(x_prev, y_prev, 0.0);
                    const geometry::point<float_t> p2(x, y, 0.0);

                    inner[0] = std::make_shared<geometry::triangle<float_t>>(center, p1, p2);
                }

                // Create outer rows
                single_glyph_t outer(2 * circle_resolution * (polygonal_resolution - 1));

                if (polygonal_resolution > 1)
                {
                    const auto circle_increment = static_cast<float_t>((2.0 * math::pi<float_t>) / circle_resolution);
                    const auto polygonal_increment = static_cast<float_t>(1.0 / polygonal_resolution);

                    auto radius_prev = polygonal_increment;

                    for (std::size_t j = 2; j <= polygonal_resolution; ++j)
                    {
                        auto radius = j * polygonal_increment;

                        auto x_prev_inner = radius_prev * std::cos(static_cast<float_t>(0.0));
                        auto y_prev_inner = radius_prev * std::sin(static_cast<float_t>(0.0));
                        auto x_prev_outer = radius * std::cos(static_cast<float_t>(0.0));
                        auto y_prev_outer = radius * std::sin(static_cast<float_t>(0.0));

                        for (std::size_t i = 1; i < circle_resolution; ++i)
                        {
                            const auto angle = i * circle_increment;

                            const auto x_inner = radius_prev * std::cos(angle);
                            const auto y_inner = radius_prev * std::sin(angle);
                            const auto x_outer = radius * std::cos(angle);
                            const auto y_outer = radius * std::sin(angle);

                            const geometry::point<float_t> p1(x_prev_inner, y_prev_inner, 0.0);
                            const geometry::point<float_t> p2(x_prev_outer, y_prev_outer, 0.0);
                            const geometry::point<float_t> p3(x_inner, y_inner, 0.0);
                            const geometry::point<float_t> p4(x_outer, y_outer, 0.0);

                            outer[(j - 2) * circle_resolution * 2 + i * 2 + 0] = std::make_shared<geometry::triangle<float_t>>(p1, p2, p4);
                            outer[(j - 2) * circle_resolution * 2 + i * 2 + 1] = std::make_shared<geometry::triangle<float_t>>(p1, p4, p3);

                            x_prev_inner = x_inner;
                            y_prev_inner = y_inner;
                            x_prev_outer = x_outer;
                            y_prev_outer = y_outer;
                        }

                        const auto x_inner = radius_prev * std::cos(static_cast<float_t>(0.0));
                        const auto y_inner = radius_prev * std::sin(static_cast<float_t>(0.0));
                        const auto x_outer = radius * std::cos(static_cast<float_t>(0.0));
                        const auto y_outer = radius * std::sin(static_cast<float_t>(0.0));

                        const geometry::point<float_t> p1(x_prev_inner, y_prev_inner, 0.0);
                        const geometry::point<float_t> p2(x_prev_outer, y_prev_outer, 0.0);
                        const geometry::point<float_t> p3(x_inner, y_inner, 0.0);
                        const geometry::point<float_t> p4(x_outer, y_outer, 0.0);

                        outer[(j - 2) * circle_resolution * 2 + 0] = std::make_shared<geometry::triangle<float_t>>(p1, p2, p4);
                        outer[(j - 2) * circle_resolution * 2 + 1] = std::make_shared<geometry::triangle<float_t>>(p1, p4, p3);

                        radius_prev = radius;
                    }
                }

                // Combine rows
                disc.insert(disc.end(), inner.begin(), inner.end());
                disc.insert(disc.end(), outer.begin(), outer.end());
            }

            // Create strips in x- and y-direction, respectively
            glyph_t strip(2 * polygonal_resolution);

            {
                const auto increment = static_cast<float_t>(1.0 / polygonal_resolution);

                auto x_prev = static_cast<float_t>(-1.0);

                const auto y_plus = strip_size / static_cast<float_t>(2.0);
                const auto y_minus = -strip_size / static_cast<float_t>(2.0);

                for (std::size_t i = 1; i <= polygonal_resolution; ++i)
                {
                    const auto x = (i * increment) * 2 - 1;

                    const geometry::point<float_t> p1(x_prev, y_plus, 0.02);
                    const geometry::point<float_t> p2(x_prev, y_minus, 0.02);
                    const geometry::point<float_t> p3(x, y_plus, 0.02);
                    const geometry::point<float_t> p4(x, y_minus, 0.02);

                    strip[(i - 1) * 2 + 0] = std::make_shared<geometry::triangle<float_t>>(p2, p4, p3);
                    strip[(i - 1) * 2 + 1] = std::make_shared<geometry::triangle<float_t>>(p2, p3, p1);

                    x_prev = x;
                }
            }

            // Return glyphs
            return std::make_tuple(disc, strip);
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

            glyph_t instance(std::get<0>(glyph_template).size() + 2 * std::get<1>(glyph_template).size());

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
                                Eigen::Matrix<float_t, 4, 1> bent_vector;
                                bent_vector << vector[0], vector[1], 
                                    scalar * (polynomial[0] * vector[0] * vector[1] + polynomial[1] * vector[0] * vector[0]
                                        + polynomial[2] * vector[1] * vector[1]) + vector[2],
                                    vector[3];

                                return bent_vector;
                            };

                            // Create object matrix
                            math::transformer<float_t, 3> trafo(translation_and_rotation * math::transformer<float_t, 3>(scale_matrix));
                            trafo.set_preprocessing(bend);

                            // Extra rotation for second strip
                            Eigen::Matrix<float_t, 4, 4> rotation_matrix;
                            rotation_matrix.setIdentity();

                            rotation_matrix(0, 0) = rotation_matrix(1, 1) = 0.0;
                            rotation_matrix(1, 0) = 1.0;
                            rotation_matrix(0, 1) = -1.0;

                            math::transformer<float_t, 3> trafo_rotation(rotation_matrix);

                            // Instantiate glyph
                            auto pos = std::transform(std::get<0>(glyph_template).begin(), std::get<0>(glyph_template).end(), instance.begin(),
                                [&trafo](const std::shared_ptr<geometry::geometric_object<float_t>> obj) { return obj->clone(trafo); });

                            pos = std::transform(std::get<1>(glyph_template).begin(), std::get<1>(glyph_template).end(), pos,
                                [&trafo](const std::shared_ptr<geometry::geometric_object<float_t>> obj) { return obj->clone(trafo); });

                            pos = std::transform(std::get<1>(glyph_template).begin(), std::get<1>(glyph_template).end(), pos,
                                [&trafo, &trafo_rotation](const std::shared_ptr<geometry::geometric_object<float_t>> obj)
                                {
                                    auto clone = obj->clone(trafo_rotation);
                                    clone->transform(trafo);
                                    return clone;
                                });

                            this->bending_glyphs->insert(instance);

                            ++num_interface_cells;
                        }
                    }
                }
            }

            // Create data array
            auto bending = std::make_shared<tpf::data::array<float_t>>("Bending");
            bending->reserve(num_interface_cells * (std::get<0>(glyph_template).size() + 2 * std::get<1>(glyph_template).size()));

            for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second && num_interface_cells != 0; ++z)
            {
                for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second && num_interface_cells != 0; ++y)
                {
                    for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second && num_interface_cells != 0; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (vof(coords) > 0 && vof(coords) < 1)
                        {
                            const auto min = bending_min(coords);
                            const auto max = bending_max(coords);

                            for (std::size_t i = 0; i < std::get<0>(glyph_template).size(); ++i)
                            {
                                bending->push_back((min + max) / 2.0);
                            }

                            for (std::size_t i = 0; i < std::get<1>(glyph_template).size(); ++i)
                            {
                                bending->push_back(min);
                            }

                            for (std::size_t i = 0; i < std::get<1>(glyph_template).size(); ++i)
                            {
                                bending->push_back(max);
                            }
                        }
                    }
                }
            }

            this->bending_glyphs->add(bending, tpf::data::topology_t::OBJECT_DATA);
        }
    }
}
