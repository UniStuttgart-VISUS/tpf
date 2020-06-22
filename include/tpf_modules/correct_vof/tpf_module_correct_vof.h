#pragma once

#include "tpf/module/tpf_module_base.h"

#include "tpf/data/tpf_grid.h"

#include <string>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for calculation of interface gradients
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class correct_vof : public module_base<
            // Input
            const data::grid<float_t, float_t, 3, 1>&,
            // Output
            data::grid<float_t, float_t, 3, 1>&,
            // Parameters
            void,
            // Callbacks
            void>
        {
        public:
            using base_t = module_base<const data::grid<float_t, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&, void, void>;

            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            correct_vof();

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
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="gradients">Output gradients</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 1>& correct_fractions) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            data::grid<float_t, float_t, 3, 1>* correct_fractions;
        };
    }
}

#include "tpf_module_correct_vof.inl"
