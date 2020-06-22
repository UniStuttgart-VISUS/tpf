#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/math/tpf_vector.h"

#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for calculation of interface curvature
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class interface_curvature : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&>,
            interface_output<data::grid<float_t, float_t, 3, 1>&>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                const data::grid<float_t, float_t, 3, 3>&>;
            using output_t = interface_output<data::grid<float_t, float_t, 3, 1>&>;

            using base_t = modules::module_base<input_t, output_t>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            interface_curvature();

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
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions,
                const data::grid<float_t, float_t, 3, 3>& gradients, const data::grid<float_t, float_t, 3, 3>& positions) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="curvature">Output Curvature</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 1>& curvature) override;

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run_algorithm() override;

        private:
            enum class direction_t
            {
                LEFT = -1,
                RIGHT = 1,
                DOWN = -2,
                UP = 2,
                FRONT = -3,
                BACK = 3
            };

            /// <summary>
            /// Calculate central first order derivative
            /// </summary>
            /// <param name="f0">Left function value</param>
            /// <param name="f1">Central function value</param>
            /// <param name="f2">Right function value</param>
            /// <param name="h0">Distance between left and center</param>
            /// <param name="h1">Distance between center and right</param>
            /// <return>First order derivative</return>
            float_t calculate_central_derivative(float_t f0, float_t f1, float_t f2, float_t h0, float_t h1) const;

            /// <summary>
            /// Calculate central second order derivative
            /// </summary>
            /// <param name="f0">Left function value</param>
            /// <param name="f1">Central function value</param>
            /// <param name="f2">Right function value</param>
            /// <param name="h0">Distance between left and center</param>
            /// <param name="h1">Distance between center and right</param>
            /// <return>Second order derivative</return>
            float_t calculate_central_second_order_derivative(float_t f0, float_t f1, float_t f2, float_t h0, float_t h1) const;

            /// <summary>
            /// Calculate the height function curvature and their respective base cells
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="direction">Most promising direction</param>
            /// <param name="fractions">VOF field</param>
            /// <return>First: curvature calculation successfull, second: [first: curvature, second: interface positions]</return>
            std::pair<bool, std::pair<float_t, std::vector<math::vec3_t<float_t>>>> height_function_curvature(const data::coords3_t& coords,
                const direction_t direction, const data::grid<float_t, float_t, 3, 1>& fractions) const;

            /// <summary>
            /// Calculate interface height for the given cell
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="direction">Most promising direction</param>
            /// <param name="fractions">VOF field</param>
            /// <return>First: height calculation successfull, second: [first: interface position, second: height]</return>
            std::pair<bool, std::pair<data::coords3_t, float_t>> calculate_interface_height(const data::coords3_t& coords,
                const direction_t direction, const data::grid<float_t, float_t, 3, 1>& fractions) const;

            /// <summary>
            /// Calculate the number of independent positions
            /// </summary>
            /// <param name="positions">Neighboring positions</param>
            /// <return>Number of independent positions</return>
            std::size_t calculate_num_independent_positions(const std::vector<math::vec3_t<float_t>>& positions, float_t delta) const;

            /// <summary>
            /// Fit a parabola through neighboring interface cells and calculate curvature
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="normal">Interface normal</param>
            /// <param name="positions">Neighboring positions</param>
            /// <param name="centroid">Own barycenter</param>
            /// <return>Curvature</return>
            float_t paraboloid_fitted_curvature(const math::vec3_t<float_t>& normal,
                const std::vector<math::vec3_t<float_t>>& positions, const math::vec3_t<float_t>& centroid) const;

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            /// Positions
            const data::grid<float_t, float_t, 3, 3>* positions;

            /// Curvature
            data::grid<float_t, float_t, 3, 1>* curvature;
        };
    }
}

#include "tpf_module_interface_curvature.inl"
