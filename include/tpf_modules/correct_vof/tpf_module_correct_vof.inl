#include "tpf_module_correct_vof.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include <string>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t correct_vof<float_t>::get_num_required_ghost_levels()
        {
            return 1;
        }

        template <typename float_t>
        inline correct_vof<float_t>::correct_vof() { }

        template <typename float_t>
        inline std::string correct_vof<float_t>::get_name() const
        {
            return std::string("Correct VOF");
        }

        template <typename float_t>
        inline void correct_vof<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions)
        {
            this->fractions = &fractions;
        }

        template <typename float_t>
        inline void correct_vof<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 1>& correct_fractions)
        {
            this->correct_fractions = &correct_fractions;
        }

        template <typename float_t>
        inline void correct_vof<float_t>::run_algorithm()
        {
            // Get input and output
            const auto& fractions = *this->fractions;

            auto& correct_fractions = *this->correct_fractions;

            // Calculate gradient at all cells
            /**
                OpenMP information

                read global:    fractions
                write global:    correct_fractions
            **/
            #pragma omp parallel for schedule(dynamic) default(none) shared(fractions, correct_fractions)
            for (auto z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                const auto z = static_cast<std::size_t>(z_omp);

                for (auto y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (auto x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        try
                        {
                            const data::coords3_t coords(x, y, z);

                            correct_fractions(coords) = fractions(coords);

                            if (fractions.is_local(coords, get_num_required_ghost_levels()) &&
                                fractions(coords) > static_cast<float_t>(0.0L) && fractions(coords) < static_cast<float_t>(1.0L))
                            {
                                // Look at neighbors
                                bool valid_neighbor = false;

                                for (long long k = -1; k <= 1; ++k)
                                {
                                    for (long long j = -1; j <= 1; ++j)
                                    {
                                        for (long long i = -1; i <= 1; ++i)
                                        {
                                            if (i != 0 || j != 0 || k != 0)
                                            {
                                                const data::coords3_t neighbor_coords = coords + data::coords3_t(i, j, k);

                                                valid_neighbor |= fractions(neighbor_coords) > static_cast<float_t>(0.0L);
                                            }
                                        }
                                    }
                                }

                                // If all neighbors are gaseous, set own fraction to zero
                                if (!valid_neighbor)
                                {
                                    correct_fractions(coords) = static_cast<float_t>(0.0L);
                                }
                            }
                        }
                        catch (const std::exception& e)
                        {
                            log::warning_message(__tpf_nested_warning_message(e.what(), "Unable to correct the VOF value for a cell."));
                        }
                        catch (...)
                        {
                            log::warning_message(__tpf_warning_message("Unable to correct the VOF value for a cell."));
                        }
                    }
                }
            }

            // Synchronize boundaries
            mpi::mpi_grid::synchronize_boundaries(correct_fractions);
        }
    }
}
