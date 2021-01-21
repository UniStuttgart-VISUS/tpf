#include "tpf_module_interface_stretching.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_linalg.h"
#include "tpf/math/tpf_matrix.h"
#include "tpf/math/tpf_vector.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include "Eigen/Dense"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t interface_stretching<float_t>::get_num_required_ghost_levels()
        {
            return 1;
        }

        template <typename float_t>
        inline interface_stretching<float_t>::interface_stretching() : fractions(nullptr), gradients(nullptr), positions(nullptr), velocities(nullptr),
            stretching(nullptr), stretching_prim_dir(nullptr), stretching_sec_dir(nullptr), stretching_absmax_dir(nullptr) { }

        template <typename float_t>
        inline std::string interface_stretching<float_t>::get_name() const
        {
            return std::string("Interface Stretching");
        }

        template <typename float_t>
        inline void interface_stretching<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
            const data::grid<float_t, float_t, 3, 3>& gradients, const data::grid<float_t, float_t, 3, 3>& positions,
            const data::grid<float_t, float_t, 3, 3>& velocities)
        {
            this->fractions = &fractions;
            this->gradients = &gradients;
            this->positions = &positions;
            this->velocities = &velocities;
        }

        template <typename float_t>
        inline void interface_stretching<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 1>& stretching,
            data::grid<float_t, float_t, 3, 3>& stretching_prim_dir, data::grid<float_t, float_t, 3, 3>& stretching_sec_dir,
            data::grid<float_t, float_t, 3, 3>& stretching_absmax_dir)
        {
            this->stretching = &stretching;
            this->stretching_prim_dir = &stretching_prim_dir;
            this->stretching_sec_dir = &stretching_sec_dir;
            this->stretching_absmax_dir = &stretching_absmax_dir;
        }

        template <typename float_t>
        inline void interface_stretching<float_t>::set_algorithm_parameters(const float_t timestep)
        {
            this->timestep = timestep;
        }

        template <typename float_t>
        inline void interface_stretching<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 3>& positions = *this->positions;
            const data::grid<float_t, float_t, 3, 3>& velocities = *this->velocities;

            data::grid<float_t, float_t, 3, 1>& stretching = *this->stretching;
            data::grid<float_t, float_t, 3, 3>& stretching_prim_dir = *this->stretching_prim_dir;
            data::grid<float_t, float_t, 3, 3>& stretching_sec_dir = *this->stretching_sec_dir;
            data::grid<float_t, float_t, 3, 3>& stretching_absmax_dir = *this->stretching_absmax_dir;

            bool good_case = true;

            // Calculate gradient at all cells
            /**
                OpenMP information

                read global:    fractions, gradients, positions, velocities
                write global:    stretching, stretching_prim_dir, stretching_sec_dir, stretching_absmax_dir
            **/
            #pragma omp parallel for schedule(dynamic) default(none) shared(fractions, gradients, positions, velocities, \
                    stretching, stretching_prim_dir, stretching_sec_dir, stretching_absmax_dir, good_case)
            for (long long z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                const auto z = static_cast<std::size_t>(z_omp);

                bool good_case_local = true;

                for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        try
                        {
                            const data::coords3_t coords(x, y, z);

                            stretching(coords) = static_cast<float_t>(0.0L);
                            stretching_prim_dir(coords).setZero();
                            stretching_sec_dir(coords).setZero();
                            stretching_absmax_dir(coords).setZero();

                            if (fractions.is_local(coords, get_num_required_ghost_levels()) &&
                                fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                            {
                                // Get position and velocity
                                std::vector<math::vec3_t<float_t>> vertices;
                                std::vector<math::vec3_t<float_t>> velocity;

                                vertices.reserve(27);
                                velocity.reserve(27);

                                vertices.push_back(positions(coords));
                                velocity.push_back(velocities(coords));

                                // Get positions and velocities at neighboring locations
                                for (long long i = -1; i <= 1; i++)
                                {
                                    for (long long j = -1; j <= 1; j++)
                                    {
                                        for (long long k = -1; k <= 1; k++)
                                        {
                                            const data::coords3_t neighbor_coords = coords + data::coords3_t(i, j, k);

                                            if (i != 0 || j != 0 || k != 0)
                                            {
                                                if (fractions(neighbor_coords) > static_cast<float_t>(0.0L))
                                                {
                                                    vertices.push_back(positions(neighbor_coords));
                                                    velocity.push_back(velocities(neighbor_coords));
                                                }
                                            }
                                        }
                                    }
                                }

                                // Calculate Jacobian matrix of gradients
                                const math::mat3_t<float_t> jacobian = math::calculate_jacobian<float_t, float_t, 3>(vertices, velocity);

                                // Calculate change in size
                                if (!jacobian.hasNaN() && !jacobian.isZero())
                                {
                                    // Calculate change in size and directions
                                    const math::vec3_t<float_t> normal = -gradients(coords).normalized();

                                    const stretching_info stretch = calculate_stretch(jacobian, normal);

                                    stretching(coords) = stretch.stretching;
                                    stretching_prim_dir(coords) = stretch.primary_direction;
                                    stretching_sec_dir(coords) = stretch.secondary_direction;
                                    stretching_absmax_dir(coords) = (std::abs(std::log2(stretch.primary_direction.norm()))
                                        > std::abs(std::log2(stretch.secondary_direction.norm())))
                                        ? stretch.primary_direction : stretch.secondary_direction;
                                }
                                else
                                {
                                    // Set dummy values
                                    const math::vec3_t<float_t> normal = -gradients(coords).normalized();
                                    const auto directions = math::orthonormal(normal);

                                    stretching(coords) = static_cast<float_t>(1.0L);
                                    stretching_prim_dir(coords) = directions.first;
                                    stretching_sec_dir(coords) = directions.second;
                                    stretching_absmax_dir(coords) = directions.first;

                                    good_case_local = false;
                                }
                            }
                        }
                        catch (const std::exception& e)
                        {
                            log::warning_message(__tpf_nested_warning_message(e.what(), "Unable to compute interface stretching for a cell."));
                        }
                        catch (...)
                        {
                            log::warning_message(__tpf_warning_message("Unable to compute interface stretching for a cell."));
                        }
                    }
                }

                #pragma omp critical(is_good_case)
                {
                    good_case &= good_case_local;
                }
            }

            if (!good_case)
            {
                log::warning_message(__tpf_warning_message("Case found with NaN- or zero-valued Jacobian."));
            }

            // Synchronize boundaries
            mpi::mpi_grid::synchronize_boundaries(stretching);
            mpi::mpi_grid::synchronize_boundaries(stretching_prim_dir);
            mpi::mpi_grid::synchronize_boundaries(stretching_sec_dir);
            mpi::mpi_grid::synchronize_boundaries(stretching_absmax_dir);
        }

        template <typename float_t>
        inline typename interface_stretching<float_t>::stretching_info interface_stretching<float_t>::calculate_stretch(const math::mat3_t<float_t>& jacobian, const math::vec3_t<float_t>& normal) const
        {
            // Multiply with time step size and add identity matrix
            math::mat3_t<float_t> scaled_jacobian = jacobian * this->timestep;
            scaled_jacobian(0, 0) += static_cast<float_t>(1.0L);
            scaled_jacobian(1, 1) += static_cast<float_t>(1.0L);
            scaled_jacobian(2, 2) += static_cast<float_t>(1.0L);

            // Create matrix consisting of two surface vectors
            const auto directions = math::orthonormal(normal);

            Eigen::Matrix<float_t, 3, 2> surface_vectors;
            surface_vectors.col(0) = directions.first;
            surface_vectors.col(1) = directions.second;

            // Create partially local jacobian
            const Eigen::Matrix<float_t, 3, 2> intermediary_jacobian = scaled_jacobian * surface_vectors;

            // Calculate local Jacobian
            const math::mat2_t<float_t> local_jacobian = intermediary_jacobian.transpose() * intermediary_jacobian;

            // Calculate Eigenvalues and Eigenvectors
            const auto eigenpair = math::calculate_eigenpair<float_t, 2>(local_jacobian);

            int prim_eigenvector = (std::abs(eigenpair.eigenvalues[0]) > std::abs(eigenpair.eigenvalues[1])) ? 0 : 1;
            int sec_eigenvector = 1 - prim_eigenvector;

            // Transform vectors back to 3D space
            math::vec3_t<float_t> prim_direction = directions.first * eigenpair.eigenvectors[prim_eigenvector][0]
                + directions.second * eigenpair.eigenvectors[prim_eigenvector][1];
            math::vec3_t<float_t> sec_direction = directions.first * eigenpair.eigenvectors[sec_eigenvector][0]
                + directions.second * eigenpair.eigenvectors[sec_eigenvector][1];

            prim_direction.normalize();
            sec_direction.normalize();

            prim_direction *= std::sqrt(eigenpair.eigenvalues[prim_eigenvector]);
            sec_direction *= std::sqrt(eigenpair.eigenvalues[sec_eigenvector]);

            // Return stretching information
            stretching_info ret;
            ret.stretching = std::sqrt(eigenpair.eigenvalues[prim_eigenvector]) * std::sqrt(eigenpair.eigenvalues[sec_eigenvector]);
            ret.primary_direction = prim_direction;
            ret.secondary_direction = sec_direction;

            return ret;
        }
    }
}
