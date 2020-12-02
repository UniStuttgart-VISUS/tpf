#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"
#include "tpf/geometry/tpf_polygon.h"
#include "tpf/geometry/tpf_polyhedron.h"

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
        class plic3 : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&>,
            interface_output<data::polydata<float_t>&>,
            interface_parameters<float_t, float_t, float_t, std::size_t, std::size_t, std::optional<float_t>>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&>;
            using output_t = interface_output<data::polydata<float_t>&>;
            using parameters_t = interface_parameters<float_t, float_t, float_t, std::size_t, std::size_t, std::optional<float_t>>;

            using base_t = module_base<input_t, output_t, parameters_t>;

        public:
            enum polygon_t {
                ERROR = 0,
                SOLID_GAS = 1,
                FLUID_GAS = 2,
                SOLID_FLUID = 3,
                SOLID_MIX = 4,
                FLUID_TYPE2 = 5,
                FLUID_TYPE13 = 6,
            };

            enum neighbor_t {
                SAME = 0,
                OTHER = 1,
                EMPTY = 2,
                MIXED = 3,
            };

            /// <summary>
            /// Constructor
            /// </summary>
            plic3();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="f">Input f-field</param>
            /// <param name="f3">Input f3-field</param>
            /// <param name="f_norm_3ph">Input normal</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& f, const data::grid<float_t, float_t, 3, 1>& f3,
                const data::grid<float_t, float_t, 3, 3>& f_norm_3ph) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="plic3_interface">Output plic3 interface</param>
            virtual void set_algorithm_output(data::polydata<float_t>& plic3_interface) override;

            /// <summary>
            /// Set parameter
            /// </summary>
            /// <param name="epsilon">Epsilon</param>
            /// <param name="error_margin_plic">Error margin for PLIC</param>
            /// <param name="error_margin_plic3">Error margin for PLIC3</param>
            /// <param name="num_iterations_plic">Number of iterations for PLIC</param>
            /// <param name="num_iterations_plic3">Number of iterations for PLIC3]</param>
            /// <param name="perturbation">Perturbation of vertices to prevent numerical instability</param>
            virtual void set_algorithm_parameters(float_t epsilon, float_t error_margin_plic, float_t error_margin_plic3,
                std::size_t num_iterations_plic, std::size_t num_iterations_plic3, std::optional<float_t> perturbation) override;

            /// <summary>Run module</summary>
            void run_algorithm() override;

        private:
            bool is_interface(const data::grid<float_t, float_t, 3, 1>& f, const data::coords3_t& coords, float_t epsilon) const;

            neighbor_t neighbor_type(const data::grid<float_t, float_t, 3, 1>& f1,
                const data::grid<float_t, float_t, 3, 1>& f2,
                const data::coords3_t& coords,
                float_t epsilon) const;

            bool has_non_zero_neighbor(const data::grid<float_t, float_t, 3, 1>& f, const data::coords3_t& coords, float_t epsilon) const;

            static std::pair<std::shared_ptr<geometry::polygon<float_t>>, float_t> reconstruct_f3_interface(const geometry::polyhedron<float_t>& poly_cell,
                float_t vof,
                const Eigen::Matrix<float_t, 3, 1>& norm,
                const Eigen::Matrix<float_t, 3, 1>& cell_coordinates,
                const Eigen::Matrix<float_t, 3, 1>& cell_size,
                float_t error_margin,
                std::size_t num_iterations,
                float_t perturbation = static_cast<float_t>(0.00001L));

            /// f-field
            const data::grid<float_t, float_t, 3, 1>* f;

            /// f3-field
            const data::grid<float_t, float_t, 3, 1>* f3;

            /// f_norm_3ph-field
            const data::grid<float_t, float_t, 3, 3>* f_norm_3ph;

            /// PLIC3 interface
            data::polydata<float_t>* plic3_interface;

            /// Epsilon
            float_t epsilon_ = static_cast<float_t>(0.00001L);

            /// Error margin plic
            float_t error_margin_plic_ = static_cast<float_t>(0.0001L);

            /// Error margin plic3
            float_t error_margin_plic3_ = static_cast<float_t>(0.0001L);

            /// Number of iterations plic
            std::size_t num_iterations_plic_ = 15;

            /// Number of iterations plic3
            std::size_t num_iterations_plic3_ = 15;

            /// Perturbation
            float_t perturbation_ = static_cast<float_t>(0.00001L);
        };
    }
}

#include "tpf_module_plic3.inl"
