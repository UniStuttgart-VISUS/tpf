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
            opt_arg<const data::grid<float_t, float_t, 3, 3>> velocities,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_min,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_max,
            opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_min,
            opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_max,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_min,
            opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_max)
        {
            this->vof = &vof;
            this->positions = &positions;

            this->velocities = get_or_default<const data::grid<float_t, float_t, 3, 3>>(velocities);

            this->stretching_min = get_or_default<const data::grid<float_t, float_t, 3, 3>>(stretching_min);
            this->stretching_max = get_or_default<const data::grid<float_t, float_t, 3, 3>>(stretching_max);

            this->bending_min = get_or_default<const data::grid<float_t, float_t, 3, 1>>(bending_min);
            this->bending_max = get_or_default<const data::grid<float_t, float_t, 3, 1>>(bending_max);
            this->bending_direction_min = get_or_default<const data::grid<float_t, float_t, 3, 3>>(bending_direction_min);
            this->bending_direction_max = get_or_default<const data::grid<float_t, float_t, 3, 3>>(bending_direction_max);
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
        inline void interface_deformation_glyph<float_t>::set_algorithm_parameters(float_t timestep)
        {
            this->timestep = timestep;
        }

        template <typename float_t>
        inline void interface_deformation_glyph<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& vof = *this->vof;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;

            // Calculate average cell size
            const auto average_cell_size = 
                std::accumulate(vof.get_cell_sizes()[0].begin(), vof.get_cell_sizes()[0].end(), 0.0) *
                std::accumulate(vof.get_cell_sizes()[0].begin(), vof.get_cell_sizes()[0].end(), 0.0) *
                std::accumulate(vof.get_cell_sizes()[0].begin(), vof.get_cell_sizes()[0].end(), 0.0);

            // Create velocity glyph
            if (this->velocities != nullptr)
            {
                const data::grid<float_t, float_t, 3, 3>& velocities = *this->velocities;

                data::polydata<float_t>& glyphs = *this->velocity_glyphs;

                // Create arrow glyph...
                const auto shaft_tip_ratio = static_cast<float_t>(0.7);
                const auto thickness_ratio = static_cast<float_t>(0.1);

                const auto shaft_thickness = thickness_ratio;
                const auto tip_thickness = static_cast<float_t>(2.0) * thickness_ratio;

                const std::size_t resolution = 16;

                const auto increment = static_cast<float_t>((2.0 * math::pi<float_t>) / resolution);

                std::vector<std::shared_ptr<geometry::geometric_object<float_t>>> arrow_glyph;

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

                // Iterate over all interface cells
                std::size_t num_interface_cells = 0;

                for (auto z = vof.get_extent()[2].first; z <= vof.get_extent()[2].second && num_interface_cells < 2; ++z)
                {
                    for (auto y = vof.get_extent()[1].first; y <= vof.get_extent()[1].second && num_interface_cells < 2; ++y)
                    {
                        for (auto x = vof.get_extent()[0].first; x <= vof.get_extent()[0].second && num_interface_cells < 2; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            if (vof(coords) > 0 && vof(coords) < 1)
                            {
                                // Get necessary information
                                const auto origin = positions(coords);
                                const auto velocity = velocities(coords);

                                const auto scale = this->timestep* velocity.norm();

                                // Create transformer
                                Eigen::Matrix<float_t, 4, 4> translation_matrix, scale_matrix, rotation_matrix;

                                translation_matrix.setIdentity();
                                translation_matrix.col(3).head(3) = origin;

                                scale_matrix.setIdentity();
                                scale_matrix *= scale;
                                scale_matrix(3, 3) = 1.0;

                                rotation_matrix.setIdentity();
                                // TODO

                                const math::transformer<float_t, 3> trafo(translation_matrix * rotation_matrix * scale_matrix);

                                // Instantiate glyph
                                auto glyph = arrow_glyph;

                                std::transform(glyph.begin(), glyph.end(), glyph.begin(),
                                    [&trafo](std::shared_ptr<geometry::geometric_object<float_t>> obj) { return obj->clone(trafo); });

                                glyphs.insert(glyph);

                                ++num_interface_cells;
                            }
                        }
                    }
                }

                // Create data array
                auto magnitudes = std::make_shared<tpf::data::array<float_t>>("Velocity Magnitude");
                magnitudes->reserve(num_interface_cells * arrow_glyph.size());

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

                                for (std::size_t i = 0; i < arrow_glyph.size(); ++i)
                                {
                                    magnitudes->push_back(magnitude);
                                }

                                --num_interface_cells; // DEBUG
                            }
                        }
                    }
                }

                glyphs.add(magnitudes, tpf::data::topology_t::OBJECT_DATA);
            }
            /*
            // Create stretching glyph
            if (this->stretching_min != nullptr && this->stretching_max != nullptr)
            {
                const data::grid<float_t, float_t, 3, 3>& stretching_min = *this->stretching_min;
                const data::grid<float_t, float_t, 3, 3>& stretching_max = *this->stretching_max;

                data::polydata<float_t>& glyphs = *this->stretching_glyphs;

                // Iterate over all interface cells
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
                                const auto stretching_1 = stretching_min(coords);
                                const auto stretching_2 = stretching_max(coords);

                                const auto scale = static_cast<float_t>(0.5 * std::pow(average_cell_size, 1.0 / 3.0));

                                // Create disc glyph

                            }
                        }
                    }
                }
            }

            // Create bending glyph
            if (this->bending_min != nullptr && this->bending_max != nullptr &&
                this->bending_direction_min != nullptr && this->bending_direction_max != nullptr)
            {
                const data::grid<float_t, float_t, 3, 1>& bending_min = *this->bending_min;
                const data::grid<float_t, float_t, 3, 1>& bending_max = *this->bending_max;
                const data::grid<float_t, float_t, 3, 3>& bending_direction_min = *this->bending_direction_min;
                const data::grid<float_t, float_t, 3, 3>& bending_direction_max = *this->bending_direction_max;

                data::polydata<float_t>& glyphs = *this->bending_glyphs;

                // Iterate over all interface cells
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
                                const auto bending_1 = bending_min(coords) * bending_direction_min(coords);
                                const auto bending_2 = bending_max(coords) * bending_direction_max(coords);

                                const auto scale = static_cast<float_t>(0.5 * std::pow(average_cell_size, 1.0 / 3.0));

                                // Create disc glyph

                            }
                        }
                    }
                }
            }*/
        }
    }
}
