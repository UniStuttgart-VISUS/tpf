#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_geometric_object.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace interface_deformation_glyph_aux
        {
            enum class arrow_size_t
            {
                dynamic, fixed_thickness, fixed_size
            };
        }

        /// <summary>
        /// Module for creating an interface deformation glyph
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class interface_deformation_glyph : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 1>>,
                opt_arg<const data::grid<float_t, float_t, 3, 1>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>>,
            interface_output<data::polydata<float_t>&, data::polydata<float_t>&, data::polydata<float_t>&>,
            interface_parameters<bool, bool, bool, float_t, interface_deformation_glyph_aux::arrow_size_t,
                float_t, float_t, int, float_t, float_t, int, int, float_t, float_t, float_t>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&,
                const data::grid<float_t, float_t, 3, 3>&,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 1>>,
                opt_arg<const data::grid<float_t, float_t, 3, 1>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>,
                opt_arg<const data::grid<float_t, float_t, 3, 3>>>;
            using output_t = interface_output<data::polydata<float_t>&, data::polydata<float_t>&, data::polydata<float_t>&>;
            using parameters_t = interface_parameters<bool, bool, bool, float_t, interface_deformation_glyph_aux::arrow_size_t,
                float_t, float_t, int, float_t, float_t, int, int, float_t, float_t, float_t>;

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
            interface_deformation_glyph();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="positions">Input volume of fluid</param>
            /// <param name="positions">Input positions</param>
            /// <param name="gradients">Input optional gradients</param>
            /// <param name="velocities">Input optional velocities</param>
            /// <param name="stretching_min">Input optional minimum stretching (direction)</param>
            /// <param name="stretching_max">Input optional maximum stretching (direction)</param>
            /// <param name="bending_min">Input optional minimum bending</param>
            /// <param name="bending_max">Input optional maximum bending</param>
            /// <param name="bending_direction_min">Input optional minimum bending direction</param>
            /// <param name="bending_direction_max">Input optional maximum bending direction</param>
            /// <param name="bending_polynomial">Input optional polynomial describing bending</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& vof,
                const data::grid<float_t, float_t, 3, 3>& positions,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> gradients,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> velocities,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_min,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> stretching_max,
                opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_min,
                opt_arg<const data::grid<float_t, float_t, 3, 1>> bending_max,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_min,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_direction_max,
                opt_arg<const data::grid<float_t, float_t, 3, 3>> bending_polynomial) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="velocity_glyphs">Output velocity glyphs</param>
            /// <param name="stretching_glyphs">Output stretching glyphs</param>
            /// <param name="bending_glyphs">Output bending glyphs</param>
            virtual void set_algorithm_output(data::polydata<float_t>& velocity_glyphs,
                data::polydata<float_t>& stretching_glyphs, data::polydata<float_t>& bending_glyphs) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="velocity_glyph">Create velocity glyphs</param>
            /// <param name="stretching_glyph">Create stretching glyphs</param>
            /// <param name="bending_glyph">Create bending glyphs</param>
            /// <param name="timestep">Time step</param>
            /// <param name="arrow_size">Sizing of arrow glyphs (velocity glyph)</param>
            /// <param name="arrow_scalar">Scalar for scaling of arrow glyphs (velocity glyph)</param>
            /// <param name="arrow_fixed_scalar">Scalar for scaling of arrow glyphs in the fixed direction (velocity glyph)</param>
            /// <param name="arrow_resolution">Resolution of arrow glyphs (velocity glyph)</param>
            /// <param name="arrow_ratio">Ratio of the shaft of arrow glyphs (velocity glyph)</param>
            /// <param name="arrow_thickness">Thickness of arrow glyphs (velocity glyph)</param>
            /// <param name="bending_disc_resolution">Resolution around the disc (bending glyph)</param>
            /// <param name="bending_polygonal_resolution">Resolution along the disc (bending glyph)</param>
            /// <param name="bending_strip_size">Thickness of the strip (bending glyph)</param>
            /// <param name="bending_size_scalar">Scalar for scaling the glyphs (bending glyph)</param>
            /// <param name="bending_scalar">Scalar for scaling the bending values (bending glyph)</param>
            virtual void set_algorithm_parameters(bool velocity_glyph, bool stretching_glyph, bool bending_glyph, float_t timestep,
                interface_deformation_glyph_aux::arrow_size_t arrow_size, float_t arrow_scalar, float_t arrow_fixed_scalar,
                int arrow_resolution, float_t arrow_ratio, float_t arrow_thickness, int bending_disc_resolution,
                int bending_polygonal_resolution, float_t bending_strip_size, float_t bending_size_scalar, float_t bending_scalar) override;

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run_algorithm() override;

        private:
            using glyph_t = std::vector<std::shared_ptr<geometry::geometric_object<float_t>>>;
            using velocity_glyph_t = glyph_t;

            using bending_glyph_t = std::tuple<glyph_t, glyph_t>;

            /// <summary>
            /// Create a velocity glyph, based at the origin and pointing in x-direction with length 1
            /// </summary>
            /// <param name="resolution">Resolution of the arrow, i.e., number of subdivisions</param>
            /// <param name="shaft_tip_ratio">Ratio of the shaft, e.g., a value 0.7 creates an arrow with shaft length 0.7 and tip length 0.3</param>
            /// <param name="thickness_ratio">Thickness, e.g., a value of 0.1 creates a shaft with radius 0.1 and a tip with radius 0.2</param>
            /// <returns>Template velocity glyph</returns>
            velocity_glyph_t create_velocity_glyph_template(std::size_t resolution, float_t shaft_tip_ratio, float_t thickness_ratio) const;

            /// <summary>
            /// "Instantiate" velocity glyphs at each interface position
            /// </summary>
            /// <param name="glyph_template">Template to instantiate</param>
            /// <param name="arrow_size">Type of scaling</param>
            /// <param name="arrow_scalar">Scaling factor</param>
            /// <param name="arrow_fixed_scalar">Scaling factor for fixed direction (thickness for fixed_thickness; all for fixed_size)</param>
            void instantiate_velocity_glyphs(const velocity_glyph_t& glyph_template,
                interface_deformation_glyph_aux::arrow_size_t arrow_size, float_t arrow_scalar, float_t arrow_fixed_scalar);

            bending_glyph_t create_bending_glyph_template(std::size_t circle_resolution, std::size_t polygonal_resolution, float_t strip_size) const;

            void instantiate_bending_glyph(const bending_glyph_t&glyph_template, float_t average_cell_size, float_t size_scalar, float_t scalar);

            /// Volume of fluid
            const data::grid<float_t, float_t, 3, 1>* vof;

            /// Positions
            const data::grid<float_t, float_t, 3, 3>* positions;

            /// Gradients
            const data::grid<float_t, float_t, 3, 3>* gradients;

            /// Velocities
            const data::grid<float_t, float_t, 3, 3>* velocities;

            /// Stretching
            const data::grid<float_t, float_t, 3, 3>* stretching_min;
            const data::grid<float_t, float_t, 3, 3>* stretching_max;

            /// Bending
            const data::grid<float_t, float_t, 3, 1>* bending_min;
            const data::grid<float_t, float_t, 3, 1>* bending_max;
            const data::grid<float_t, float_t, 3, 3>* bending_direction_min;
            const data::grid<float_t, float_t, 3, 3>* bending_direction_max;
            const data::grid<float_t, float_t, 3, 3>* bending_polynomial;

            /// Glyphs
            data::polydata<float_t>* velocity_glyphs;
            data::polydata<float_t>* stretching_glyphs;
            data::polydata<float_t>* bending_glyphs;

            /// Create glyphs?
            bool velocity_glyph;
            bool stretching_glyph;
            bool bending_glyph;

            /// Time step for scaling the velocities
            float_t timestep;

            /// Properties of the velocity glyphs
            interface_deformation_glyph_aux::arrow_size_t arrow_size;
            float_t arrow_scalar;
            float_t arrow_fixed_scalar;
            int arrow_resolution;
            float_t arrow_ratio;
            float_t arrow_thickness;

            /// Properties of the bending glyph
            int bending_disc_resolution;
            int bending_polygonal_resolution;
            float_t bending_strip_size;
            float_t bending_size_scalar;
            float_t bending_scalar;

            /// Properties of the stretching glyph
            
        };
    }
}

#include "tpf_module_interface_deformation_glyph.inl"
