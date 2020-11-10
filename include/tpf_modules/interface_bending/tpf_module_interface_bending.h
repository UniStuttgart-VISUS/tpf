#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/algorithm/tpf_least_squares.h"

#include "tpf/data/tpf_grid.h"

#include "tpf/math/tpf_vector.h"

#include <string>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for calculation of interface bending
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class interface_bending : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 3>&>,
            interface_output<data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
                data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&>,
            interface_parameters<float_t>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 3>&>;
            using output_t = interface_output<data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
                data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&>;
            using parameters_t = interface_parameters<float_t>;

            using base_t = module_base<input_t, output_t, parameters_t>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            interface_bending();
            
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
            /// <param name="positions">Input positions</param>
            /// <param name="velocities">Input velocities</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions, const data::grid<float_t, float_t, 3, 3>& gradients,
                const data::grid<float_t, float_t, 3, 3>& positions, const data::grid<float_t, float_t, 3, 3>& velocities) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="min_curvature">Output minimum curvature change</param>
            /// <param name="max_curvature">Output maximum curvature change</param>
            /// <param name="absmax_curvature">Output absolute maximum curvature change</param>
            /// <param name="min_curvature_dir">Output min direction of curvature change</param>
            /// <param name="max_curvature_dir">Output max direction of curvature change</param>
            /// <param name="absmax_curvature_dir">Output absolute max direction of curvature change</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 1>& min_curvature, data::grid<float_t, float_t, 3, 1>& max_curvature,
                data::grid<float_t, float_t, 3, 1>& absmax_curvature, data::grid<float_t, float_t, 3, 3>& min_curvature_dir,
                data::grid<float_t, float_t, 3, 3>& max_curvature_dir, data::grid<float_t, float_t, 3, 3>& absmax_curvature_dir) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="timestep">Time step delta</param>
            virtual void set_algorithm_parameters(float_t timestep) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Struct for storing curvature information
            /// </summary>
            struct curvature_info
            {
                float_t first_curvature;
                float_t second_curvature;

                math::vec2_t<float_t> first_direction;
                math::vec2_t<float_t> second_direction;
            };

            /// <summary>
            /// Calculate the principal curvatures and their directions using the fundamental forms
            /// of a polynomial of the form a_xy xy + a_xx x^2 + a_yy y^2
            /// </summary>
            /// <param name="poly">Polynomial</param>
            /// <return>Curvatures and fundamental forms</return>
            curvature_info calculate_curvature(const algorithm::polynomial<float_t>& poly) const;

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            // Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            // Barycenters
            const data::grid<float_t, float_t, 3, 3>* positions;

            // Velocities
            const data::grid<float_t, float_t, 3, 3>* velocities;

            /// Minimum curvature change
            data::grid<float_t, float_t, 3, 1>* min_curvature;

            /// Maximum curvature change
            data::grid<float_t, float_t, 3, 1>* max_curvature;

            /// Absolute maximum curvature change
            data::grid<float_t, float_t, 3, 1>* absmax_curvature;

            /// Min direction of curvature change
            data::grid<float_t, float_t, 3, 3>* min_curvature_dir;

            /// Max direction of curvature change
            data::grid<float_t, float_t, 3, 3>* max_curvature_dir;

            /// Absolute max direction of curvature change
            data::grid<float_t, float_t, 3, 3>* absmax_curvature_dir;

            // Time step delta
            float_t timestep;
        };
    }
}

#include "tpf_module_interface_bending.inl"
