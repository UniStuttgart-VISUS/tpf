#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"

#include "tpf/math/tpf_matrix.h"
#include "tpf/math/tpf_vector.h"

#include <string>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for calculation of interface stretching
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class interface_stretching : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 3>&>,
            interface_output<data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 3>&,
                data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&>,
            interface_parameters<float_t>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&, const data::grid<float_t, float_t, 3, 3>&>;
            using output_t = interface_output<data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 3>&,
                data::grid<float_t, float_t, 3, 3>&, data::grid<float_t, float_t, 3, 3>&>;
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
            interface_stretching();

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
            /// <param name="stretching">Output stretching</param>
            /// <param name="stretching_prim_dir">Output primary direction of stretching</param>
            /// <param name="stretching_sec_dir">Output secondary direction of stretching</param>
            /// <param name="stretching_absmax_dir">Output absolute max direction of stretching</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 1>& stretching, data::grid<float_t, float_t, 3, 3>& stretching_prim_dir,
                data::grid<float_t, float_t, 3, 3>& stretching_sec_dir, data::grid<float_t, float_t, 3, 3>& stretching_absmax_dir) override;

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
            /// Struct for storing stretching information
            /// </summary>
            struct stretching_info
            {
                float_t stretching;
                math::vec3_t<float_t> primary_direction;
                math::vec3_t<float_t> secondary_direction;
            };

            /// <summary>
            /// Calculate stretch using a two-dimensional metric tensor
            /// </summary>
            /// <param name="jacobian">Jacobian matrix</param>
            /// <param name="normal">Normal</param>
            /// <return>[stretching factor, [primary direction of stretching, secondary direction of stretching]]</return>
            stretching_info calculate_stretch(const math::mat3_t<float_t>& jacobian, const math::vec3_t<float_t>& normal) const;

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            // Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            // Barycenters
            const data::grid<float_t, float_t, 3, 3>* positions;

            // Velocities
            const data::grid<float_t, float_t, 3, 3>* velocities;

            /// Stretching
            data::grid<float_t, float_t, 3, 1>* stretching;

            // Primary direction of stretching
            data::grid<float_t, float_t, 3, 3>* stretching_prim_dir;

            // Secondary direction of stretching
            data::grid<float_t, float_t, 3, 3>* stretching_sec_dir;

            // Absolute max direction of stretching
            data::grid<float_t, float_t, 3, 3>* stretching_absmax_dir;

            // Time step delta
            float_t timestep;
        };
    }
}

#include "tpf_module_interface_stretching.inl"
