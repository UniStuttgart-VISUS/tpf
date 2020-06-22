#pragma once

#include "tpf/module/tpf_module_base.h"

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
            // Input
            std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 1>&,
            const data::grid<float_t, float_t, 3, 3>&>,
            // Output
            data::polydata<float_t>&,
            // Parameters
            std::tuple<std::optional<std::size_t>, std::optional<std::size_t>>,
            // Callbacks
            void>
        {
            using base_t = module_base<std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&>, data::polydata<float_t>&, std::tuple<std::optional<std::size_t>, std::optional<std::size_t>>, void>;

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
            /// <param name="input">Input [f-field, f3-field, f_norm_3ph]</param>
            virtual void set_algorithm_input(std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&> input) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="plic3_interface">Output plic3 interface</param>
            virtual void set_algorithm_output(data::polydata<float_t>& plic3_interface) override;

            /// <summary>
            /// Set parameter
            /// </summary>
            /// <param name="parameters">[Number of iterations for PLIC, Number of iterations for PLIC3]</param>
            virtual void set_algorithm_parameters(std::tuple<std::optional<std::size_t>, std::optional<std::size_t>> parameters) override;

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

            std::size_t plic_iterations_ = 15;
            std::size_t plic3_iterations_ = 15;
        };
    }
}

#include "tpf_module_plic3.inl"
