#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "Eigen/Dense"

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
        class interface_gradient : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&>,
            interface_output<data::grid<float_t, float_t, 3, 3>&>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&>;
            using output_t = interface_output<data::grid<float_t, float_t, 3, 3>&>;

            using base_t = module_base<input_t, output_t>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            interface_gradient();

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
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 3>& gradients) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            template<typename T>
            friend class plic3;

            /// <summary>
            /// Calculate the gradient at the given position
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="fractions">VOF</param>
            /// <return>Gradient</return>
            static Eigen::Matrix<float_t, 3, 1> calculate_gradient(const data::coords3_t& coords, const data::grid<float_t, float_t, 3, 1>& fractions);

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            data::grid<float_t, float_t, 3, 3>* gradients;
        };
    }
}

#include "tpf_module_interface_gradient.inl"
