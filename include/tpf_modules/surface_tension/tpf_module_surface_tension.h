#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

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
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&>,
            interface_output<data::grid<float_t, float_t, 3, 3>&>,
            interface_parameters<float_t, float_t, float_t>>
        {
        public:
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 1>&>;
            using output_t = interface_output<data::grid<float_t, float_t, 3, 3>&>;
            using parameters_t = interface_parameters<float_t, float_t, float_t>;

            using base_t = module_base<input_t, output_t, parameters_t>;

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
            /// <param name="fractions">Input fractions</param>
            /// <param name="gradients">Input gradients</param>
            /// <param name="curvature">Input curvature</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
                const data::grid<float_t, float_t, 3, 3>& gradients, const data::grid<float_t, float_t, 3, 1>& curvature) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="surface_tension_force">Output surface tension</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 3>& surface_tension_force) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="coefficient">Coefficient depending on the media involved</param>
            /// <param name="density">Density</param>
            /// <param name="timestep">Time step</param>
            virtual void set_algorithm_parameters(float_t coefficient, float_t density, float_t timestep) override;

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
