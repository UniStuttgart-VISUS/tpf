#pragma once

#include "tpf/module/tpf_module_base.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_plane.h"

#include "tpf/math/tpf_vector.h"

#include <functional>
#include <optional>
#include <string>
#include <tuple>
#include <utility>

namespace tpf
{
    namespace modules
    {
        namespace fluid_position_aux
        {
            /// <summary>
            /// Position type
            /// </summary>
            enum class position_t
            {
                INVALID = -1,
                CELL_CENTER, FLUID_CENTER, INTERFACE
            };
        }

        /// <summary>
        /// Module for calculation of interface positions
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class fluid_position : public module_base<
            // Input
            std::tuple<const data::grid<float_t, float_t, 3, 1>&, std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>,
            // Output
            std::tuple<data::grid<float_t, float_t, 3, 3>&, data::polydata<float_t>&>,
            // Parameters
            fluid_position_aux::position_t,
            // Callbacks
            void>
        {
            using base_t = module_base<std::tuple<const data::grid<float_t, float_t, 3, 1>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>,
                std::tuple<data::grid<float_t, float_t, 3, 3>&, data::polydata<float_t>&>, fluid_position_aux::position_t, void>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            fluid_position();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="input">Input [fractions, optional gradients]</param>
            virtual void set_algorithm_input(std::tuple<const data::grid<float_t, float_t, 3, 1>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>> input) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="output">Output [positions grid, positions points]</param>
            virtual void set_algorithm_output(std::tuple<data::grid<float_t, float_t, 3, 3>&, data::polydata<float_t>&> output) override;

            /// <summary>
            /// Set paramter
            /// </summary>
            /// <param name="position_type">Position type</param>
            virtual void set_algorithm_parameters(fluid_position_aux::position_t position_type) override;

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run_algorithm() override;

        private:
            /// <summary>
            /// Calculate the barycenter of the interface in the given cell
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="fraction">VOF</param>
            /// <param name="normal">Interface normal</param>
            /// <param name="cell_sizes">Cell sizes for every direction</param>
            /// <return>Barycenter</return>
            math::vec3_t<float_t> calculate_interface_barycenter(float_t fraction, const math::vec3_t<float_t>& normal,
                const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const;

            /// <summary>
            /// Calculate the barycenter of the fluid in the given cell
            /// </summary>
            /// <param name="coords">Coordinates</param>
            /// <param name="fraction">VOF</param>
            /// <param name="normal">Interface normal</param>
            /// <param name="cell_sizes">Cell sizes for every direction</param>
            /// <return>Barycenter</return>
            math::vec3_t<float_t> calculate_fluid_barycenter(float_t fraction, const math::vec3_t<float_t>& normal,
                const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const;

            /// <summary>
            /// Calculate PLIC interface
            /// </summary>
            /// <param name="fraction">Cell VOF</param>
            /// <param name="normal">Cell normal</param>
            /// <param name="cell_sizes">Cell size</param>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <returns>PLIC and corner</returns>
            std::pair<geometry::plane<float_t>, math::vec3_t<float_t>> calculate_plic(float_t fraction, const math::vec3_t<float_t>& normal,
                const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const;

            /// <summary>
            /// Calculate cell boundaries
            /// </summary>
            /// <param name="cell_sizes">Cell size</param>
            /// <param name="cell_coordinates">Cell coordinates</param>
            /// <returns>Cell cuboid</returns>
            geometry::cuboid<float_t> calculate_cell(const math::vec3_t<float_t>& cell_sizes, const math::vec3_t<float_t>& cell_coordinates) const;

            /// <summary>
            /// Calculate the distance between the innermost corner to the interface
            /// </summary>
            /// <author>Grzegorz K. Karch</author>
            /// <param name="fraction">VOF</param>
            /// <param name="normal">Interface normal</param>
            /// <param name="cell_sizes">Cell sizes for every direction</param>
            /// <return>Distance</return>
            float_t calculate_distance_to_interface(float_t fraction, const math::vec3_t<float_t>& normal, const math::vec3_t<float_t>& cell_sizes) const;

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            /// Positions
            data::grid<float_t, float_t, 3, 3>* positions_grid;
            data::polydata<float_t>* positions_points;

            /// Position type
            fluid_position_aux::position_t position_type;
        };
    }
}

#include "tpf_module_fluid_position.inl"
