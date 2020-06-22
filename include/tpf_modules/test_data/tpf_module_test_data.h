#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"

#include <functional>
#include <optional>
#include <string>

namespace tpf
{
    namespace modules
    {
        namespace test_data_aux
        {
            enum geometry_t
            {
                INVALID = -1,
                ORB, RING
            };
        }

        /// <summary>
        /// Module for creation of test data
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class test_data : public module_base<
            interface_output<data::grid<float_t, float_t, 3, 1>&, std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>>,
            interface_parameters<test_data_aux::geometry_t, std::size_t, bool, bool, bool, bool,
                std::optional<float_t>, std::optional<float_t>, std::optional<float_t>, std::optional<float_t>>>
        {
        public:
            using output_t = interface_output<data::grid<float_t, float_t, 3, 1>&, std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>>;
            using parameters_t = interface_parameters<test_data_aux::geometry_t, std::size_t, bool, bool, bool, bool,
                std::optional<float_t>, std::optional<float_t>, std::optional<float_t>, std::optional<float_t>>;

            using base_t = module_base<output_t, parameters_t>;

            /// <summary>
            /// Constructor
            /// </summary>
            test_data();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="fractions">Output volume of fluid field</param>
            /// <param name="velocities">Optional output velocity field</param>
            virtual void set_algorithm_output(data::grid<float_t, float_t, 3, 1>& fractions,
                std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>> velocities) override;

            /// <summary>
            /// Set paramters
            /// </summary>
            /// <param name="test_geometry">Test geometry_t</param>
            /// <param name="num_cells">Number of cells per direction</param>
            /// <param name="spinning">Add spinning velocities?</param>
            /// <param name="moving">Add moving velocities?</param>
            /// <param name="expanding">Add expanding velocities?</param>
            /// <param name="rotating">Add rotating velocities?</param>
            /// <param name="spinning_magnitude">Magnitude of spinning velocity</param>
            /// <param name="moving_magnitude">Magnitude of movement velocity</param>
            /// <param name="expanding_magnitude">Magnitude of expansion velocity</param>
            /// <param name="rotating_magnitude">Magnitude of rotation velocity</param>
            virtual void set_algorithm_parameters(test_data_aux::geometry_t test_geometry, std::size_t num_cells, bool spinning,
                bool moving, bool expanding, bool rotating, std::optional<float_t> spinning_magnitude, std::optional<float_t> moving_magnitude,
                std::optional<float_t> expanding_magnitude, std::optional<float_t> rotating_magnitude ) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Create a square cartesian grid
            /// </summary>
            /// <param name="num_cells">Number of cells per direction</param>
            /// <returns>Grids</returns>
            void create_square_cartesian_grids(const std::size_t num_cells);

            /// Fractions
            data::grid<float_t, float_t, 3, 1>* fractions;

            /// Velocities
            data::grid<float_t, float_t, 3, 3>* velocities;

            /// Test geometry_t
            test_data_aux::geometry_t test_geometry;

            /// Number of cells per direction
            std::size_t num_cells;

            /// Velocities
            bool spinning;
            bool moving;
            bool expanding;
            bool rotating;

            /// Velocity magnitudes
            float_t spinning_magnitude;
            float_t moving_magnitude;
            float_t expanding_magnitude;
            float_t rotating_magnitude;
        };
    }
}

#include "tpf_module_test_data.inl"
