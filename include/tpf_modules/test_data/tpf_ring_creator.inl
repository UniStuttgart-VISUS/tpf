#include "tpf_ring_creator.h"

#include "tpf_geometry_creator.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_halfspace.h"
#include "tpf/geometry/tpf_intersection.h"
#include "tpf/geometry/tpf_plane.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_polyhedron.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_vector.h"

#include <cstddef>

namespace tpf
{
    namespace modules
    {
        namespace test_data_aux
        {
            template <typename floatp_t>
            inline ring_creator<floatp_t>::ring_creator(data::grid<floatp_t, floatp_t, 3, 1>* fractions,
                data::grid<floatp_t, floatp_t, 3, 3>* velocities, std::size_t num_cells)
                : geometry_creator<floatp_t>(fractions, velocities, num_cells),
                ring_center(static_cast<floatp_t>(1.0L / 2.0L), static_cast<floatp_t>(1.0L / 2.0L), static_cast<floatp_t>(1.0L / 2.0L)),
                ring_plane_normal(math::vec3_t<floatp_t>(static_cast<floatp_t>(1.0L), static_cast<floatp_t>(1.0L), static_cast<floatp_t>(0.0L)).normalized()),
                ring_inner_radius(static_cast<floatp_t>(1.0L / 5.0L)), ring_outer_radius(static_cast<floatp_t>(1.0L / 4.0L)),
                ring_radius((ring_outer_radius - ring_inner_radius) / static_cast<floatp_t>(2.0L))
            {
#ifdef __tpf_detailed
                log::info_message(__tpf_info_message("Creating torus with inner radius ", ring_inner_radius, " and outer radius ", ring_outer_radius));
#endif

                create();
            }

            template <typename floatp_t>
            inline void ring_creator<floatp_t>::add_spinning(const floatp_t magnitude)
            {
                const auto& fractions = *this->fractions;
                auto& velocities = *this->velocities;

                for (std::size_t z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
                {
                    for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                    {
                        for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            // Set spinning velocities
                            if (fractions(coords) > static_cast<floatp_t>(0.0L))
                            {
                                const math::vec3_t<floatp_t> cell_center = fractions.get_cell_coordinates(coords);

                                velocities(coords) += this->ring_plane_normal.cross(cell_center - this->ring_center) * magnitude;
                            }
                        }
                    }
                }
            }

            template <typename floatp_t>
            inline void ring_creator<floatp_t>::add_movement(const floatp_t magnitude)
            {
                const auto& fractions = *this->fractions;
                auto& velocities = *this->velocities;

                // Set movement direction
                const math::vec3_t<floatp_t> movement_direction(static_cast<floatp_t>(1.0L), static_cast<floatp_t>(1.0L), static_cast<floatp_t>(1.0L));

                for (std::size_t z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
                {
                    for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                    {
                        for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            // Set movement velocities
                            if (fractions(coords) > static_cast<floatp_t>(0.0L))
                            {
                                velocities(coords) += movement_direction.normalized() * magnitude;
                            }
                        }
                    }
                }
            }

            template <typename floatp_t>
            inline void ring_creator<floatp_t>::add_expansion(const floatp_t magnitude)
            {
                const auto& fractions = *this->fractions;
                auto& velocities = *this->velocities;

                for (std::size_t z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
                {
                    for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                    {
                        for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            // Set expansion velocities
                            if (fractions(coords) > static_cast<floatp_t>(0.0L))
                            {
                                const math::vec3_t<floatp_t> cell_center = fractions.get_cell_coordinates(coords);
                                const math::vec3_t<floatp_t> normal = (cell_center - this->ring_center).normalized();

                                velocities(coords) += normal.normalized() * magnitude;
                            }
                        }
                    }
                }
            }

            template <typename floatp_t>
            inline void ring_creator<floatp_t>::add_rotation(const floatp_t magnitude)
            {
                const auto& fractions = *this->fractions;
                auto& velocities = *this->velocities;

                for (std::size_t z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
                {
                    for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                    {
                        for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            // Set rotation velocities
                            if (fractions(coords) > static_cast<floatp_t>(0.0L))
                            {
                                const math::vec3_t<floatp_t> cell_center = fractions.get_cell_coordinates(coords);
                                const math::vec3_t<floatp_t> ring_to_cell_center = (cell_center - this->ring_center).normalized();
                                const floatp_t distance = (cell_center - this->ring_center).norm();

                                const math::vec3_t<floatp_t> nearest_ring_center = this->ring_plane_normal.cross(ring_to_cell_center).cross(this->ring_plane_normal).normalized()
                                    * ((this->ring_outer_radius + this->ring_inner_radius) / static_cast<floatp_t>(2.0L)) + this->ring_center;

                                const math::vec3_t<floatp_t> plic_normal = (cell_center - nearest_ring_center).normalized();

                                velocities(coords) += ring_to_cell_center.cross(this->ring_plane_normal).cross(plic_normal).normalized() * magnitude;
                            }
                        }
                    }
                }
            }

            template <typename floatp_t>
            inline void ring_creator<floatp_t>::create()
            {
                auto& fractions = *this->fractions;
                auto& velocities = *this->velocities;

                for (std::size_t z = fractions.get_extent()[2].first; z <= fractions.get_extent()[2].second; ++z)
                {
                    for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                    {
                        for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                        {
                            const data::coords3_t coords(x, y, z);

                            this->positions(coords).setZero();

                            // Create VOF field
                            const math::vec3_t<floatp_t> cell_center = fractions.get_cell_coordinates(coords);
                            const math::vec3_t<floatp_t> ring_to_cell_center = (cell_center - this->ring_center).normalized();
                            const floatp_t distance = (cell_center - this->ring_center).norm();
                            const floatp_t diagonal = fractions.get_cell_sizes(coords).norm();

                            if (distance - diagonal > this->ring_outer_radius || distance + diagonal < this->ring_inner_radius)
                            {
                                fractions(coords) = static_cast<floatp_t>(0.0L);
                            }
                            else
                            {
                                const math::vec3_t<floatp_t> node_coords = fractions.get_node_coordinates(coords);
                                const math::vec3_t<floatp_t> cell_size = fractions.get_cell_sizes(coords);

                                const math::vec3_t<floatp_t> nearest_ring_center = this->ring_plane_normal.cross(ring_to_cell_center).cross(this->ring_plane_normal).normalized()
                                    * ((this->ring_outer_radius + this->ring_inner_radius) / static_cast<floatp_t>(2.0L)) + this->ring_center;

                                const math::vec3_t<floatp_t> plic_normal = (cell_center - nearest_ring_center).normalized();

                                const geometry::plane<floatp_t> plic(geometry::point<floatp_t>(nearest_ring_center + plic_normal * this->ring_radius), plic_normal);
                                const geometry::cuboid<floatp_t> box(geometry::point<floatp_t>(node_coords), geometry::point<floatp_t>(node_coords + cell_size));

                                auto intersections = geometry::intersect_with(plic, box);

                                if (intersections.size() < 3)
                                {
                                    if (geometry::is_in_positive_halfspace(geometry::point<floatp_t>(cell_center), plic))
                                    {
                                        fractions(coords) = static_cast<floatp_t>(0.0L);
                                    }
                                    else
                                    {
                                        fractions(coords) = static_cast<floatp_t>(1.0L);
                                        this->positions(coords) = cell_center;
                                    }
                                }
                                else
                                {
                                    for (const auto& corner : box.get_points())
                                    {
                                        if (!geometry::is_in_positive_halfspace(geometry::point<floatp_t>(corner), plic))
                                        {
                                            intersections.push_back(geometry::point<floatp_t>(corner));
                                        }
                                    }

                                    const geometry::polyhedron<floatp_t> polyhedron(intersections);

                                    const floatp_t cuboid_volume = box.calculate_volume();
                                    const floatp_t polyhedron_volume = polyhedron.calculate_volume();

                                    fractions(coords) = polyhedron_volume / cuboid_volume;
                                    this->positions(coords) = polyhedron.calculate_centroid().get_points()[0];
                                }
                            }

                            // Create velocities
                            velocities(coords).setZero();
                        }
                    }
                }
            }
        }
    }
}
