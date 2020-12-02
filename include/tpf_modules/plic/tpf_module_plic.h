#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_polygon.h"

#include "Eigen/Dense"

#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for interface reconstruction using PLIC
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class plic : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&>,
            interface_output<data::polydata<float_t>&>,
            interface_parameters<float_t, std::size_t, std::optional<float_t>>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&>;
            using output_t = interface_output<data::polydata<float_t>&>;
            using parameters_t = interface_parameters<float_t, std::size_t, std::optional<float_t>>;

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
            plic();

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
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
                const data::grid<float_t, float_t, 3, 3>& gradients) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="plic_interface">Output plic interface</param>
            virtual void set_algorithm_output(data::polydata<float_t>& plic_interface) override;

            /// <summary>
            /// Set parameter
            /// </summary>
            /// <param name="error_margin">Error margin</param>
            /// <param name="num_iterations">Number of iterations</param>
            /// <param name="perturbation">Perturbation of vertices to prevent numerical instability</param>
            virtual void set_algorithm_parameters(float_t error_margin, std::size_t num_iterations, std::optional<float_t> perturbation) override;

            /// <summary>Run module</summary>
            virtual void run_algorithm() override;

        private:
            template<typename T>
            friend class plic3;

            /// <summary>
            /// Compute ISO value for a given VOF value
            /// </summary>
            /// <param name="vof">Volume of fluid</param>
            /// <param name="gradient">Gradient</param>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <param name="cell_size">Cell size</param>
            /// <param name="num_iterations">Number of iterations</param>
            /// <param name="perturbation">Perturbation of vertices to prevent numerical instability</param>
            /// <returns>PLIC interface</returns>
            static std::pair<std::shared_ptr<geometry::polygon<float_t>>, float_t> reconstruct_interface(float_t vof,
                const Eigen::Matrix<float_t, 3, 1>& gradient, const Eigen::Matrix<float_t, 3, 1>& cell_coordinates,
                const Eigen::Matrix<float_t, 3, 1>& cell_size, float_t error_margin, std::size_t num_iterations,
                float_t perturbation = static_cast<float_t>(0.00001L));

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            /// PLIC interface
            data::polydata<float_t>* plic_interface;

            /// Error margin
            float_t error_margin;

            /// Number of iterations
            std::size_t num_iterations;

            /// Perturbation
            float_t perturbation;
        };
    }
}

#include "tpf_module_plic.inl"
