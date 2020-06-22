#include "tpf_module_surface_tension.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include <string>
#include <tuple>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t surface_tension<float_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t>
        inline surface_tension<float_t>::surface_tension() { }

        template <typename float_t>
        inline std::string surface_tension<float_t>::get_name() const
        {
            return std::string("Surface Tension");
        }

        template <typename float_t>
        inline void surface_tension<float_t>::set_algorithm_input(std::tuple<const data::grid<float_t, float_t, 3, 1>&,
            const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&> input)
        {
            this->fractions = &std::get<0>(input);
            this->gradients = &std::get<1>(input);
            this->curvature = &std::get<2>(input);
        }

        template <typename float_t>
        inline void surface_tension<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 3>& surface_tension_force)
        {
            this->surface_tension_force = &surface_tension_force;
        }

        template <typename float_t>
        inline void surface_tension<float_t>::set_algorithm_parameters(std::tuple<float_t, float_t, float_t> parameters)
        {
            this->coefficient = std::get<0>(parameters);
            this->density = std::get<1>(parameters);
            this->timestep = std::get<2>(parameters);
        }

        template <typename float_t>
        inline void surface_tension<float_t>::run_algorithm()
        {
            // Get input and output
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;
            const data::grid<float_t, float_t, 3, 3>& gradients = *this->gradients;
            const data::grid<float_t, float_t, 3, 1>& curvature = *this->curvature;

            data::grid<float_t, float_t, 3, 3>& surface_tension_force = *this->surface_tension_force;

            // Calculate gradient at all cells
            /**
                OpenMP information

                read global:    fractions, gradients, curvature
                write global:    surface_tension_force
            **/
            #pragma omp parallel for schedule(dynamic) default(none) shared(fractions, gradients, curvature, surface_tension_force)
            for (long long z_omp = static_cast<long long>(fractions.get_extent()[2].first); z_omp <= static_cast<long long>(fractions.get_extent()[2].second); ++z_omp)
            {
                const std::size_t z = static_cast<std::size_t>(z_omp);

                for (std::size_t y = fractions.get_extent()[1].first; y <= fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = fractions.get_extent()[0].first; x <= fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        surface_tension_force(coords).setZero();

                        if (fractions.is_local(coords, get_num_required_ghost_levels()) &&
                            fractions(coords) > static_cast<float_t>(0.0) && fractions(coords) < static_cast<float_t>(1.0))
                        {
                            surface_tension_force(coords) = gradients(coords) * curvature(coords) * this->coefficient * this->density * this->timestep;
                        }
                    }
                }
            }
        }
    }
}
