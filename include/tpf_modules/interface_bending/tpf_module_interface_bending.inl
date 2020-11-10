#include "tpf_module_interface_bending.h"

#include "tpf/algorithm/tpf_least_squares.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_transformer.h"
#include "tpf/math/tpf_vector.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t interface_bending<float_t>::get_num_required_ghost_levels()
        {
            return 2;
        }

        template <typename float_t>
        inline interface_bending<float_t>::interface_bending() : fractions(nullptr), gradients(nullptr), positions(nullptr), velocities(nullptr),
            min_curvature(nullptr), max_curvature(nullptr), absmax_curvature(nullptr), min_curvature_dir(nullptr), max_curvature_dir(nullptr),
            absmax_curvature_dir(nullptr) { }

        template <typename float_t>
        inline std::string interface_bending<float_t>::get_name() const
        {
            return std::string("Interface Bending");
        }

        template <typename float_t>
        inline void interface_bending<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
            const data::grid<float_t, float_t, 3, 3>& gradients, const data::grid<float_t, float_t, 3, 3>& positions,
            const data::grid<float_t, float_t, 3, 3>& velocities)
        {
            this->fractions = &fractions;
            this->gradients = &gradients;
            this->positions = &positions;
            this->velocities = &velocities;
        }

        template <typename float_t>
        inline void interface_bending<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 1>& min_curvature,
            data::grid<float_t, float_t, 3, 1>& max_curvature, data::grid<float_t, float_t, 3, 1>& absmax_curvature,
            data::grid<float_t, float_t, 3, 3>& min_curvature_dir, data::grid<float_t, float_t, 3, 3>& max_curvature_dir,
            data::grid<float_t, float_t, 3, 3>& absmax_curvature_dir)
        {
            this->min_curvature = &min_curvature;
            this->max_curvature = &max_curvature;
            this->absmax_curvature = &absmax_curvature;
            this->min_curvature_dir = &min_curvature_dir;
            this->max_curvature_dir = &max_curvature_dir;
            this->absmax_curvature_dir = &absmax_curvature_dir;
        }

        template <typename float_t>
        inline void interface_bending<float_t>::set_algorithm_parameters(const float_t timestep)
        {
            this->timestep = timestep;
        }

        template <typename float_t>
        inline void interface_bending<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& velocities = *this->velocities;

            data::grid<float_t, float_t, 3, 1>& min_curvature = *this->min_curvature;
            data::grid<float_t, float_t, 3, 1>& max_curvature = *this->max_curvature;
            data::grid<float_t, float_t, 3, 1>& absmax_curvature = *this->absmax_curvature;
            data::grid<float_t, float_t, 3, 3>& min_curvature_dir = *this->min_curvature_dir;
            data::grid<float_t, float_t, 3, 3>& max_curvature_dir = *this->max_curvature_dir;
            data::grid<float_t, float_t, 3, 3>& absmax_curvature_dir = *this->absmax_curvature_dir;

            // Calculate bending at all cells
            /**
                OpenMP information

                read global:    fractions, gradients, positions, velocities
                write global:    min_curvature, max_curvature, absmax_curvature, min_curvature_dir, max_curvature_dir, absmax_curvature_dir
            **/
            #pragma omp parallel for schedule(dynamic) default(none) shared(fractions, gradients, positions, velocities, \
                    min_curvature, max_curvature, absmax_curvature, min_curvature_dir, max_curvature_dir, absmax_curvature_dir)
            for (long long z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                std::size_t z = static_cast<std::size_t>(z_omp);

                for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        min_curvature(coords) = static_cast<float_t>(0.0L);
                        max_curvature(coords) = static_cast<float_t>(0.0L);
                        absmax_curvature(coords) = static_cast<float_t>(0.0L);
                        min_curvature_dir(coords).setZero();
                        max_curvature_dir(coords).setZero();
                        absmax_curvature_dir(coords).setZero();

                        if (fractions.is_local(coords, get_num_required_ghost_levels()) &&
                            fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                        {
                            // Get neighbouring barycenters and advect them
                            std::vector<math::vec3_t<float_t>> points, velocity;
                            points.reserve(27);
                            velocity.reserve(27);

                            for (long long k = -2; k <= 2; ++k)
                            {
                                for (long long j = -2; j <= 2; ++j)
                                {
                                    for (long long i = -2; i <= 2; ++i)
                                    {
                                        const data::coords3_t neighbour_pos = coords + data::coords3_t(i, j, k);

                                        if (fractions(neighbour_pos) > static_cast<float_t>(0.0L) && fractions(neighbour_pos) < static_cast<float_t>(1.0L))
                                        {
                                            points.push_back(positions(neighbour_pos));
                                            velocity.push_back(velocities(neighbour_pos));
                                        }
                                    }
                                }
                            }

                            // Get central point and its normal
                            const math::vec3_t<float_t> normal = -gradients(coords).normalized();

                            const math::vec3_t<float_t> origin = positions(coords);
                            const math::vec3_t<float_t> origin_velocity = velocities(coords);

                            // Subtract velocity of central point from all others
                            for (auto& vel : velocity)
                            {
                                vel -= origin_velocity;
                            }

                            // Calculate angular velocity
                            math::vec3_t<float_t> angular_velocity;
                            angular_velocity.setZero();

                            for (std::size_t i = 0; i < points.size(); ++i)
                            {
                                const math::vec3_t<float_t> dist = points[i] - origin;

                                if (!dist.isZero())
                                {
                                    angular_velocity += dist.cross(velocity[i]) / dist.squaredNorm();
                                }
                            }

                            angular_velocity /= static_cast<float_t>(points.size() - 1);

                            // Eliminate angular velocity
                            for (std::size_t i = 0; i < points.size(); ++i)
                            {
                                const math::vec3_t<float_t> dist = points[i] - origin;

                                if (!dist.isZero())
                                {
                                    velocity[i] -= angular_velocity.cross(dist);
                                }
                            }

                            // Advect with simulation velocity or surface tension
                            std::vector<math::vec3_t<float_t>> advected_points;
                            advected_points.reserve(points.size());

                            for (std::size_t i = 0; i < points.size(); ++i)
                            {
                                advected_points.push_back(points[i] + this->timestep * velocity[i]);
                            }

                            const math::vec3_t<float_t> advected_origin = origin + this->timestep * origin_velocity;

                            // Transform points to 2D
                            math::transformer<float_t, 3> trafo(origin, normal);
                            math::transformer<float_t, 3> trafo_advected(advected_origin, normal);

                            for (std::size_t i = 0; i < points.size(); ++i)
                            {
                                trafo.transform_inverse_inplace(points[i]);
                                trafo_advected.transform_inverse_inplace(advected_points[i]);
                            }

                            // Calculate paraboloids
                            const auto paraboloid = algorithm::least_squares(points.begin(), points.end());
                            const auto paraboloid_advected = algorithm::least_squares(advected_points.begin(), advected_points.end());

                            const auto difference = paraboloid_advected - paraboloid;

                            // Calculate curvature of difference
                            auto curvature_difference = calculate_curvature(difference);

                            // Sort curvatures and corresponding direcions
                            if (curvature_difference.first_curvature > curvature_difference.second_curvature)
                            {
                                auto temp_curv = curvature_difference.first_curvature;
                                auto temp_curv_dir = curvature_difference.first_direction;

                                curvature_difference.first_curvature = curvature_difference.second_curvature;
                                curvature_difference.first_direction = curvature_difference.second_direction;

                                curvature_difference.second_curvature = temp_curv;
                                curvature_difference.second_direction = temp_curv_dir;
                            }

                            // Transform back to 3D
                            const math::transformer<float_t, 3, 0> trafo_vec(math::vec3_t<float_t>(), normal);

                            math::vec3_t<float_t> curv_dir_1, curv_dir_2;
                            curv_dir_1 << curvature_difference.first_direction, static_cast<float_t>(0.0L);
                            curv_dir_2 << curvature_difference.second_direction, static_cast<float_t>(0.0L);

                            trafo_vec.transform_inplace(curv_dir_1);
                            trafo_vec.transform_inplace(curv_dir_2);

                            // Set output
                            min_curvature(coords) = curvature_difference.first_curvature;
                            max_curvature(coords) = curvature_difference.second_curvature;
                            absmax_curvature(coords) = (std::abs(min_curvature(coords)) > std::abs(max_curvature(coords))) ? min_curvature(coords) : max_curvature(coords);

                            min_curvature_dir(coords) = curv_dir_1;
                            max_curvature_dir(coords) = curv_dir_2;
                            absmax_curvature_dir(coords) = (std::abs(min_curvature(coords)) > std::abs(max_curvature(coords))) ? curv_dir_1 : curv_dir_2;
                        }
                    }
                }
            }

            // Synchronize boundaries
            mpi::mpi_grid::synchronize_boundaries(min_curvature);
            mpi::mpi_grid::synchronize_boundaries(max_curvature);
            mpi::mpi_grid::synchronize_boundaries(absmax_curvature);
            mpi::mpi_grid::synchronize_boundaries(min_curvature_dir);
            mpi::mpi_grid::synchronize_boundaries(max_curvature_dir);
            mpi::mpi_grid::synchronize_boundaries(absmax_curvature_dir);
        }

        template <typename float_t>
        inline typename interface_bending<float_t>::curvature_info interface_bending<float_t>::calculate_curvature(const algorithm::polynomial<float_t>& poly) const
        {
            // Calculate mean and Gaussian curvature
            const float_t meanCurv = poly.a_xx + poly.a_yy;
            const float_t gaussCurv = static_cast<float_t>(4.0L) * poly.a_xx * poly.a_yy - poly.a_xy * poly.a_xy;

            // Calculate principal curvature and direction
            curvature_info curvatures;

            if (math::equals(poly.a_xy, static_cast<float_t>(0.0L)))
            {
                curvatures.first_curvature = static_cast<float_t>(2.0L) * poly.a_xx;
                curvatures.second_curvature = static_cast<float_t>(2.0L) * poly.a_yy;

                curvatures.first_direction << static_cast<float_t>(1.0L), static_cast<float_t>(0.0L);
                curvatures.second_direction << static_cast<float_t>(0.0L), static_cast<float_t>(1.0L);
            }
            else
            {
                curvatures.first_curvature = meanCurv - std::sqrt(meanCurv * meanCurv - gaussCurv);
                curvatures.second_curvature = meanCurv + std::sqrt(meanCurv * meanCurv - gaussCurv);

                curvatures.first_direction << poly.a_xy, curvatures.first_curvature - static_cast<float_t>(2.0L) * poly.a_xx;
                curvatures.second_direction << poly.a_xy, curvatures.second_curvature - static_cast<float_t>(2.0L) * poly.a_xx;
            }

            return curvatures;
        }
    }
}
