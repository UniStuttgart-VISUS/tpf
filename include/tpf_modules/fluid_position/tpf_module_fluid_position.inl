#include "tpf_module_fluid_position.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_geometric_object.h"
#include "tpf/geometry/tpf_intersection.h"
#include "tpf/geometry/tpf_plane.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_polygon.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_transformer.h"
#include "tpf/math/tpf_vector.h"

#include <algorithm>
#include <cmath>
#include <functional>
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
        inline std::size_t fluid_position<float_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t>
        inline fluid_position<float_t>::fluid_position() { }

        template <typename float_t>
        inline std::string fluid_position<float_t>::get_name() const
        {
            return std::string("Fluid Position");
        }

        template <typename float_t>
        inline void fluid_position<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
            std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>> gradients)
        {
            this->fractions = &fractions;
            this->gradients = get_or_default<const data::grid<float_t, float_t, 3, 3>>(gradients);
        }

        template <typename float_t>
        inline void fluid_position<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 3>& positions_grid, data::polydata<float_t>& positions_points)
        {
            this->positions_grid = &positions_grid;
            this->positions_points = &positions_points;
        }

        template <typename float_t>
        inline void fluid_position<float_t>::set_algorithm_parameters(const fluid_position_aux::position_t position_type)
        {
            this->position_type = position_type;
        }

        template <typename float_t>
        inline void fluid_position<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;

            data::grid<float_t, float_t, 3, 3>& positions_grid = *this->positions_grid;
            data::polydata<float_t>& positions_points = *this->positions_points;

            std::vector<std::shared_ptr<geometry::geometric_object<float_t>>> positions_points_thread;

            /**
                OpenMP information

                write critical:    positions_points
                write private:     positions_points_thread

                read global:       fractions
                write global:      positions_grid
            **/
            #pragma omp parallel for schedule(dynamic) default(none) private(positions_points_thread) shared(fractions, positions_grid, positions_points)
            for (long long z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                const std::size_t z = static_cast<std::size_t>(z_omp);

                for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        try
                        {
                            const data::coords3_t coords(x, y, z);

                            positions_grid(coords).setZero();

                            if (fractions.is_local(coords, get_num_required_ghost_levels()) && fractions(coords) > static_cast<float_t>(0.0L))
                            {
                                if (this->position_type == fluid_position_aux::position_t::CELL_CENTER)
                                {
                                    positions_grid(coords) = fractions.get_cell_coordinates(coords);
                                    positions_points_thread.push_back(std::make_shared<geometry::point<float_t>>(positions_grid(coords)));
                                }
                                else
                                {
                                    const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
                                    const math::vec3_t<float_t> normal = -gradients(coords).normalized();

                                    if (this->position_type == fluid_position_aux::position_t::FLUID_CENTER)
                                    {
                                        if (fractions(coords) < static_cast<float_t>(1.0L))
                                        {
                                            positions_grid(coords) = calculate_fluid_barycenter(fractions(coords), normal,
                                                fractions.get_cell_sizes(coords), fractions.get_cell_coordinates(coords));

                                            positions_points_thread.push_back(std::make_shared<geometry::point<float_t>>(positions_grid(coords)));
                                        }
                                        else
                                        {
                                            positions_grid(coords) = fractions.get_cell_coordinates(coords);
                                            positions_points_thread.push_back(std::make_shared<geometry::point<float_t>>(positions_grid(coords)));
                                        }
                                    }
                                    else if (this->position_type == fluid_position_aux::position_t::INTERFACE)
                                    {
                                        if (fractions(coords) < static_cast<float_t>(1.0L))
                                        {
                                            positions_grid(coords) = calculate_interface_barycenter(fractions(coords), normal,
                                                fractions.get_cell_sizes(coords), fractions.get_cell_coordinates(coords));

                                            positions_points_thread.push_back(std::make_shared<geometry::point<float_t>>(positions_grid(coords)));
                                        }
                                    }
                                }
                            }
                        }
                        catch (const std::exception& e)
                        {
                            log::warning_message(__tpf_nested_warning_message(e.what(), "Unable to compute fluid position for a cell."));
                        }
                        catch (...)
                        {
                            log::warning_message(__tpf_warning_message("Unable to compute fluid position for a cell."));
                        }
                    }
                }

                // Merge private with global positions
                #pragma omp critical(write_positions)
                {
                    for (auto positions : positions_points_thread)
                    {
                        positions_points.insert(positions);
                    }
                }

                positions_points_thread.clear();
            }
        }

        template <typename float_t>
        inline math::vec3_t<float_t> fluid_position<float_t>::calculate_interface_barycenter(const float_t fraction, const math::vec3_t<float_t>& normal,
            const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const
        {
            // Define plic plane, cuboid representing the cell and intersect them
            const std::pair<geometry::plane<float_t>, math::vec3_t<float_t>> plic = calculate_plic(fraction, normal, cell_sizes, cell_coordinates);
            const geometry::cuboid<float_t> cell = calculate_cell(cell_sizes, cell_coordinates);

            std::vector<geometry::point<float_t>> intersection_points = geometry::template intersect_with<float_t>(plic.first, cell);

            // Return the corner as interface position when there are no intersections
            if (intersection_points.size() < 3)
            {
                return plic.second;
            }

            // Calculate centroid
            return geometry::polygon<float_t>(intersection_points, normal, true).calculate_centroid();
        }

        template <typename float_t>
        inline math::vec3_t<float_t> fluid_position<float_t>::calculate_fluid_barycenter(const float_t fraction, const math::vec3_t<float_t>& normal,
            const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const
        {
            // Define plic plane, cuboid representing the cell and intersect them
            const std::pair<geometry::plane<float_t>, math::vec3_t<float_t>> plic = calculate_plic(fraction, normal, cell_sizes, cell_coordinates);
            const geometry::cuboid<float_t> cell = calculate_cell(cell_sizes, cell_coordinates);

            std::vector<geometry::point<float_t>> intersection_points = geometry::template intersect_with<float_t>(plic.first, cell);

            // Return the corner as interface position when there are no intersections
            if (intersection_points.size() < 4)
            {
                return calculate_interface_barycenter(fraction, normal, cell_sizes, cell_coordinates);
            }

            // Create polyhedron and calculate its centroid
            try
            {
                geometry::polyhedron<float_t> fluid(intersection_points);

                return fluid.calculate_centroid().get_vertex();
            }
#ifdef __tpf_debug
            catch (const std::runtime_error& e)
#else
            catch (const std::runtime_error&)
#endif
            {
#ifdef __tpf_debug
                log::warning_message(__tpf_nested_warning_message(e.what(), "Unable to calculate fluid barycenter. Using interface barycenter instead."));
#endif

                return calculate_interface_barycenter(fraction, normal, cell_sizes, cell_coordinates);
            }
        }

        template <typename float_t>
        inline std::pair<geometry::plane<float_t>, math::vec3_t<float_t>> fluid_position<float_t>::calculate_plic(const float_t fraction, const math::vec3_t<float_t>& normal,
            const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const
        {
            // Calculate distance between interface and corner on the inside
            const auto distance = calculate_distance_to_interface(fraction, normal, cell_sizes);

            const math::vec3_t<float_t> corner(
                cell_coordinates[0] + ((normal[0] > static_cast<float_t>(0.0L)) ? static_cast<float_t>(-1.0L) : static_cast<float_t>(1.0L)) * cell_sizes[0] / static_cast<float_t>(2.0L),
                cell_coordinates[1] + ((normal[1] > static_cast<float_t>(0.0L)) ? static_cast<float_t>(-1.0L) : static_cast<float_t>(1.0L)) * cell_sizes[1] / static_cast<float_t>(2.0L),
                cell_coordinates[2] + ((normal[2] > static_cast<float_t>(0.0L)) ? static_cast<float_t>(-1.0L) : static_cast<float_t>(1.0L)) * cell_sizes[2] / static_cast<float_t>(2.0L));

            // Define plic plane, cuboid representing the cell and intersect them
            const math::vec3_t<float_t> point = corner + distance * normal;
            const geometry::plane<float_t> plic(geometry::point<float_t>(point), normal);

            return std::make_pair(plic, corner);
        }

        template <typename float_t>
        inline geometry::cuboid<float_t> fluid_position<float_t>::calculate_cell(const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const
        {
            const geometry::point<float_t> min(cell_coordinates - cell_sizes / static_cast<float_t>(2.0L));
            const geometry::point<float_t> max(cell_coordinates + cell_sizes / static_cast<float_t>(2.0L));

            const geometry::cuboid<float_t> cell(min, max);

            return cell;
        }

        template <typename float_t>
        inline float_t fluid_position<float_t>::calculate_distance_to_interface(const float_t fraction, const math::vec3_t<float_t>& normal,
            const math::vec3_t<float_t>& cell_sizes) const
        {
            const std::size_t MAXITER = 100;

            std::size_t ih, i1, i2, i3;
            float_t d1, d2, d3;
            float_t n1, n2, n3;
            float_t nd[3], nd1, nd2, nd3, ndsum;
            float_t volume;
            float_t li, lii, liii, liv, lv;
            float_t la, lb, ll, sla, slb, sumla, sumlb, dlstar;
            std::size_t niter;
            float_t d2d3, n2rd3, n2n3;
            float_t lstar;
            float_t epsf = static_cast<float_t>(0.001L);

            niter = 0;
            nd[0] = std::abs(normal[0] * cell_sizes[0]);
            nd[1] = std::abs(normal[1] * cell_sizes[1]);
            nd[2] = std::abs(normal[2] * cell_sizes[2]);
            i1 = 0;
            i2 = 1;
            i3 = 2;

            // [i3] < [i2] < [i1] (?)
            if (nd[0] < nd[1])
            {
                i1 = 1;
                i2 = 0;
            }
            if (nd[i2] < nd[2])
            {
                i3 = i2;
                i2 = 2;
            }
            if (nd[i1] < nd[i2])
            {
                ih = i1;
                i1 = i2;
                i2 = ih;
            }

            d1 = cell_sizes[i1];
            d2 = cell_sizes[i2];
            d3 = cell_sizes[i3];

            n1 = std::abs(normal[i1]);
            n2 = std::abs(normal[i2]);
            n3 = std::abs(normal[i3]);

            nd1 = nd[i1];
            nd2 = nd[i2];
            nd3 = nd[i3];

            ndsum = nd1 + nd2 + nd3;
            volume = d1 * d2 * d3;

            d2d3 = d2 * d3;
            n2rd3 = n2 / d3;
            n2n3 = n2 * n3;

            if (fraction < static_cast<float_t>(0.000001L))
                return 0.0;
            if (fraction > static_cast<float_t>(0.999999L))
                return ndsum;

            li = nd1;
            lii = nd1 + nd3;
            liii = nd1 + nd2;
            liv = liii + nd3;
            lv = liv + nd1;

            dlstar = static_cast<float_t>(0.0L);
            lstar = static_cast<float_t>(0.5L) * liv;

            sumla = static_cast<float_t>(0.0L);
            sumlb = static_cast<float_t>(0.5L) * volume * n1;

            while (1)
            {
                if (fabs((sumlb - sumla) / volume / n1 - fraction) < epsf || niter > MAXITER)
                    break;
                ++niter;
                la = lstar;
                lb = lstar + nd1;

                //-------------------------------------------------------------------------
                // calculation of ll, sla, slb
                // Area 1
                if (la >= static_cast<float_t>(0.0L) && la <= li)
                {
                    sla = static_cast<float_t>(0.0L);
                    sumla = static_cast<float_t>(0.0L);
                }
                else if (lb >= static_cast<float_t>(0.0L) && lb <= li)
                {
                    slb = static_cast<float_t>(0.0L);
                    sumlb = static_cast<float_t>(0.0L);
                }
                // Area 2
                if (la >= li && la <= lii)
                {
                    ll = la - li;
                    sla = ll * ll / (static_cast<float_t>(2.0L)*n2n3);
                    sumla = ll * ll*ll / (static_cast<float_t>(6.0L)*n2n3);
                }
                else if (lb >= li && lb <= lii)
                {
                    ll = lb - li;
                    slb = ll * ll / (static_cast<float_t>(2.0L)*n2n3);
                    sumlb = ll * ll*ll / (static_cast<float_t>(6.0L)*n2n3);
                }
                // Area 3
                if (la >= lii && la <= liii)
                {
                    ll = la - lii;
                    sla = nd3 / n2rd3 / static_cast<float_t>(2.0L) + ll / n2rd3;
                    sumla = (static_cast<float_t>(3.0L)*ll*(nd3 + ll) + nd3 * nd3) / n2rd3 / static_cast<float_t>(6.0L);
                }
                else if (lb >= lii && lb <= liii)
                {
                    ll = lb - lii;
                    slb = nd3 / n2rd3 / static_cast<float_t>(2.0L) + ll / n2rd3;
                    sumlb = (static_cast<float_t>(3.0L)*ll*(nd3 + ll) + nd3 * nd3) / n2rd3 / static_cast<float_t>(6.0L);
                }
                // Area 4
                if (la >= liii && la <= liv)
                {
                    ll = liv - la;
                    sla = d2d3 - static_cast<float_t>(0.5L)*ll*ll / n2n3;
                    sumla = (static_cast<float_t>(3.0L)*nd2*nd3*nd3 + static_cast<float_t>(3.0L)*nd2*nd2*nd3 -
                        static_cast<float_t>(6.0L)*nd2*nd3*ll + ll * ll*ll) / n2n3 / static_cast<float_t>(6.0L);
                }
                else if (lb >= liii && lb <= liv)
                {
                    ll = liv - lb;
                    slb = d2d3 - static_cast<float_t>(0.5L)*ll*ll / n2n3;
                    sumlb = (static_cast<float_t>(3.0L)*nd2*nd3*nd3 + static_cast<float_t>(3.0L)*nd2*nd2*nd3 -
                        static_cast<float_t>(6.0L)*nd2*nd3*ll + ll * ll*ll) / n2n3 / static_cast<float_t>(6.0L);
                }
                // Area 5
                if (la >= liv && la <= lv)
                {
                    ll = la - liv;
                    sla = d2d3;
                    sumla = d2d3 * (ll + static_cast<float_t>(0.5L)*nd2 + static_cast<float_t>(0.5L)*nd3);
                }
                else if (lb >= liv && lb <= lv)
                {
                    ll = lb - liv;
                    slb = d2d3;
                    sumlb = d2d3 * (ll + static_cast<float_t>(0.5L)*nd2 + static_cast<float_t>(0.5L)*nd3);
                }
                //-------------------------------------------------------------------------
                //      dlstar = (sumlb - sumla - f*volume*n1) / MAX((slb - sla), emf*d2d3);
                dlstar = (sumlb - sumla - fraction * volume*n1) / (slb - sla);
                lstar = lstar - dlstar;
                lstar = std::max(lstar, static_cast<float_t>(0.0L));
                lstar = std::min(lstar, liv);
            }

            return lstar;
        }
    }
}
