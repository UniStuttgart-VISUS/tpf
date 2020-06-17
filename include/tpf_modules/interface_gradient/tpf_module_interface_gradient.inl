#include "tpf_module_interface_gradient.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_stencil.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include "Eigen/Dense"

#include <stdexcept>
#include <string>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t interface_gradient<float_t>::get_num_required_ghost_levels()
        {
            return 1;
        }

        template <typename float_t>
        inline interface_gradient<float_t>::interface_gradient() : fractions(nullptr), gradients(nullptr) { }

        template <typename float_t>
        inline std::string interface_gradient<float_t>::get_name() const
        {
            return std::string("Interface gradient");
        }

        template <typename float_t>
        inline void interface_gradient<float_t>::set_input(const data::grid<float_t, float_t, 3, 1>& fractions)
        {
            this->fractions = &fractions;
        }

        template <typename float_t>
        inline void interface_gradient<float_t>::set_output(data::grid<float_t, float_t, 3, 3>& gradients)
        {
            this->gradients = &gradients;
        }

        template <typename float_t>
        inline void interface_gradient<float_t>::run_algorithm()
        {
            if (this->fractions == nullptr || this->gradients == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Module run started before assigning input and output."));
            }

            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;

            data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;

            // Calculate gradient at all cells
            /**
                OpenMP information

                read global:    fractions
                write global:    gradients
            **/
            #pragma omp parallel for schedule(dynamic) default(none) shared(fractions, gradients)
            for (auto z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                const auto z = static_cast<std::size_t>(z_omp);

                // Calculate gradient
                for (auto y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (auto x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        gradients(coords).setZero();

                        if (fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                        {
                            gradients(coords) = calculate_gradient(coords, fractions);
                        }
                    }
                }
            }

            // Synchronize boundaries
            mpi::mpi_grid::synchronize_boundaries(gradients);
        }

        template <typename float_t>
        inline Eigen::Matrix<float_t, 3, 1> interface_gradient<float_t>::calculate_gradient(const data::coords3_t& coords, const data::grid<float_t, float_t, 3, 1>& fractions)
        {
            // Create 3x3 stencil around the coordinates
            using stencil_t = data::stencil<const float_t, float_t, 3, 1>;
            const stencil_t stencil(fractions, coords, 3); // stencil type does not matter here, we test for is_on_grid before reading

            // Compute gradient on nodes
            Eigen::Matrix<float_t, 3, 1> gradient;
            gradient.setZero();

            float_t number_of_points = static_cast<float_t>(0.0L);

            float_t dfm1, dfm2, dc_sum_inv_2;

            for (std::size_t k = 0; k < 2; ++k)
            {
                const std::size_t km = k;
                const std::size_t kp = k + 1;

                for (std::size_t j = 0; j < 2; ++j) {
                    const std::size_t jm = j;
                    const std::size_t jp = j + 1;

                    for (std::size_t i = 0; i < 2; ++i) {
                        const std::size_t im = i;
                        const std::size_t ip = i + 1;

                        const data::coords3_t smaller(im, jm, km);
                        const data::coords3_t larger(ip, jp, kp);

                        const auto smaller_cell_size = stencil.get_cell_sizes(smaller);
                        const auto larger_cell_size = stencil.get_cell_sizes(larger);

                        const Eigen::Matrix<float_t, 3, 1> dc = (smaller_cell_size + larger_cell_size);
                        dc_sum_inv_2 = static_cast<float_t>(2L) / (dc[0] * dc[1] * dc[2]);

                        const data::coords3_t c[8] = {
                            data::coords3_t(im, jm, km),
                            data::coords3_t(ip, jm, km),
                            data::coords3_t(im, jp, km),
                            data::coords3_t(ip, jp, km),
                            data::coords3_t(im, jm, kp),
                            data::coords3_t(ip, jm, kp),
                            data::coords3_t(im, jp, kp),
                            data::coords3_t(ip, jp, kp)
                        };

                        const bool on_grid[8] = {
                            stencil.is_on_grid(c[0]),
                            stencil.is_on_grid(c[1]),
                            stencil.is_on_grid(c[2]),
                            stencil.is_on_grid(c[3]),
                            stencil.is_on_grid(c[4]),
                            stencil.is_on_grid(c[5]),
                            stencil.is_on_grid(c[6]),
                            stencil.is_on_grid(c[7])
                        };

                        if (!on_grid[0] || !on_grid[1] || !on_grid[2] || !on_grid[3] || !on_grid[4] || !on_grid[5] || !on_grid[6] || !on_grid[7]) {
                            continue;
                        }

                        number_of_points += static_cast<float_t>(1.0L);

                        const float_t f[8] = {
                            stencil(c[0]),
                            stencil(c[1]),
                            stencil(c[2]),
                            stencil(c[3]),
                            stencil(c[4]),
                            stencil(c[5]),
                            stencil(c[6]),
                            stencil(c[7])
                        };

                        dfm1 = (f[7] - f[6]) * smaller_cell_size[2] + (f[3] - f[2]) * larger_cell_size[2];
                        dfm2 = (f[5] - f[4]) * smaller_cell_size[2] + (f[1] - f[0]) * larger_cell_size[2];
                        const float_t nx = (dfm1 * smaller_cell_size[1] + dfm2 * larger_cell_size[1]) * dc_sum_inv_2;

                        dfm1 = (f[7] - f[5]) * smaller_cell_size[2] + (f[3] - f[1]) * larger_cell_size[2];
                        dfm2 = (f[6] - f[4]) * smaller_cell_size[2] + (f[2] - f[0]) * larger_cell_size[2];
                        const float_t ny = (dfm1 * smaller_cell_size[0] + dfm2 * larger_cell_size[0]) * dc_sum_inv_2;

                        dfm1 = (f[7] - f[3]) * smaller_cell_size[1] + (f[5] - f[1]) * larger_cell_size[1];
                        dfm2 = (f[6] - f[2]) * smaller_cell_size[1] + (f[4] - f[0]) * larger_cell_size[1];
                        const float_t nz = (dfm1 * smaller_cell_size[0] + dfm2 * larger_cell_size[0]) * dc_sum_inv_2;

                        // Gradient points from f inwards
                        gradient += Eigen::Matrix<float_t, 3, 1>(nx, ny, nz);
                    }
                }
            }

            if (number_of_points > static_cast<float_t>(0.0L)) {
                return gradient / number_of_points;
            } else {
                return gradient;
            }
        }
    }
}