#pragma once

#include "tpf/module/tpf_module_base.h"

#include "tpf/data/tpf_grid.h"

#include <string>
#include <tuple>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for calculation of surface tension
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class surface_tension : public module_base<
            // Input
            std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&>,
            // Output
            data::grid<float_t, float_t, 3, 3>&,
            // Parameters
            std::tuple<float_t, float_t, float_t>,
            // Callbacks
            void>
        {
        public:
            using base_t = module_base<std::tuple<const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&>,
                data::grid<float_t, float_t, 3, 3>&, std::tuple<float_t, float_t, float_t>, void>;

            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            surface_tension();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="input">Input [fractions, gradients, curvature]</param>
            virtual void set_algorithm_input(std::tuple<const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&> input) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="surface_tension_force">Output surface tension</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 3>& surface_tension_force) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="parameters">Parameters: [coefficient depending on the media involved, density, time step]</param>
            virtual void set_algorithm_parameters(std::tuple<float_t, float_t, float_t> parameters) override;

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run_algorithm() override;

        private:
            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            /// Curvature
            const data::grid<float_t, float_t, 3, 1>* curvature;

            /// Surface tension
            data::grid<float_t, float_t, 3, 3>* surface_tension_force;

            /// Coefficient depending on the media involved
            float_t coefficient;

            /// Density
            float_t density;

            /// Time step
            float_t timestep;
        };
    }
}

#include "tpf_module_surface_tension.inl"
