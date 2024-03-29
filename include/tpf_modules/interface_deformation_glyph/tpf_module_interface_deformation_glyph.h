#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_geometric_object.h"
#include "tpf/geometry/tpf_mesh.h"

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
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

            template <typename float_t>
            struct velocity_params_t
            {
                /// Size type
                arrow_size_t size;

                /// Scalar for scaling of the glyphs
                float_t scalar;

                /// Scalar for scaling of the glyphs in the fixed direction
                float_t fixed_scalar;

                /// Resolution
                int resolution;

                /// Length of the shaft
                float_t shaft_length;

                /// Thickness of the arrow
                float_t shaft_thickness;
            };

            template <typename float_t>
            struct stretching_params_t
            {
                /// Scalar for increasing the stretching effect
                float_t scalar;

                /// Sharpness of superellipsoid
                float_t sharpness;

                /// Scalar for scaling the glyphs
                float_t size_scalar;

                /// Radius of the whole within the disc
                float_t hole_radius;

                /// Offset in normal direction
                float_t offset;

                /// Resolution around the disc
                int disc_resolution;

                /// Disc bending for pseudo-3D effect
                float_t disc_bending;

                /// Thickness of the strip
                bool show_strip;
                float_t strip_size;

                /// Thickness of the reference circle
                bool show_reference;
                float_t reference_size;

                /// z-offset of the strips and the reference circle to prevent z-fighting
                float_t z_offset;
            };

            template <typename float_t>
            struct bending_params_t
            {
                /// Scalar for scaling the bending values
                float_t scalar;

                /// Scalar for scaling the glyphs
                float_t size_scalar;

                /// Offset in normal direction
                float_t offset;

                /// Resolution around the disc
                int disc_resolution;

                /// Resolution along the disc
                int polynomial_resolution;

                /// Thickness of the strip
                bool show_strip;
                float_t strip_size;

                /// z-offset of the strips to prevent z-fighting
                float_t z_offset;
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
            interface_output<data::polydata<float_t>&,
                std::array<std::reference_wrapper<data::polydata<float_t>>, 3>,
                std::array<std::reference_wrapper<data::polydata<float_t>>, 2>>,
            interface_parameters<bool, bool, bool, float_t,
                interface_deformation_glyph_aux::velocity_params_t<float_t>,
                interface_deformation_glyph_aux::stretching_params_t<float_t>,
                interface_deformation_glyph_aux::bending_params_t<float_t>>>
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
            using output_t = interface_output<data::polydata<float_t>&,
                std::array<std::reference_wrapper<data::polydata<float_t>>, 3>,
                std::array<std::reference_wrapper<data::polydata<float_t>>, 2>>;
            using parameters_t = interface_parameters<bool, bool, bool, float_t,
                interface_deformation_glyph_aux::velocity_params_t<float_t>,
                interface_deformation_glyph_aux::stretching_params_t<float_t>,
                interface_deformation_glyph_aux::bending_params_t<float_t>>;

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
                std::array<std::reference_wrapper<data::polydata<float_t>>, 3> stretching_glyphs,
                std::array<std::reference_wrapper<data::polydata<float_t>>, 2> bending_glyphs) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="velocity_glyph">Create velocity glyphs</param>
            /// <param name="stretching_glyph">Create stretching glyphs</param>
            /// <param name="bending_glyph">Create bending glyphs</param>
            /// <param name="timestep">Time step</param>
            /// <param name="velocity_parameters">Parameters for the velocity glyph</param>
            /// <param name="stretching_parameters">Parameters for the stretching glyph</param>
            /// <param name="bending_parameters">Parameters for the bending glyph</param>
            virtual void set_algorithm_parameters(bool velocity_glyph, bool stretching_glyph, bool bending_glyph, float_t timestep,
                interface_deformation_glyph_aux::velocity_params_t<float_t> velocity_parameters,
                interface_deformation_glyph_aux::stretching_params_t<float_t> stretching_parameters,
                interface_deformation_glyph_aux::bending_params_t<float_t> bending_parameters) override;

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run_algorithm() override;

        private:
            using glyph_t = std::shared_ptr<geometry::mesh<float_t>>;
            using glyph_tex_t = std::pair<glyph_t, std::shared_ptr<std::vector<Eigen::Matrix<float_t, 2, 1>>>>;
            using velocity_glyph_t = glyph_t;
            using stretching_glyph_t = std::tuple<glyph_tex_t, glyph_tex_t, glyph_tex_t, glyph_tex_t>;
            using bending_glyph_t = std::tuple<glyph_t, glyph_t, glyph_t>;

            /// <summary>
            /// Create a velocity glyph, based at the origin and pointing in x-direction with length 1
            /// </summary>
            /// <param name="resolution">Resolution of the arrow, i.e., number of subdivisions</param>
            /// <param name="shaft_length">Ratio of the shaft, e.g., a value 0.7 creates an arrow with shaft length 0.7 and tip length 0.3</param>
            /// <param name="shaft_thickness">Thickness, e.g., a value of 0.1 creates a shaft with radius 0.1 and a tip with radius 0.2</param>
            /// <returns>Template velocity glyph</returns>
            velocity_glyph_t create_velocity_glyph_template(std::size_t resolution, float_t shaft_length, float_t shaft_thickness) const;

            /// <summary>
            /// "Instantiate" velocity glyphs at each interface position
            /// </summary>
            /// <param name="glyph_template">Template to instantiate</param>
            /// <param name="average_cell_size">Average cell size of the data</param>
            /// <param name="arrow_size">Type of scaling</param>
            /// <param name="arrow_scalar">Scaling factor</param>
            /// <param name="arrow_fixed_scalar">Scaling factor for fixed direction (thickness for fixed_thickness; all for fixed_size)</param>
            void instantiate_velocity_glyphs(const velocity_glyph_t& glyph_template, float_t average_cell_size,
                interface_deformation_glyph_aux::arrow_size_t arrow_size, float_t arrow_scalar, float_t arrow_fixed_scalar);

            /// <summary>
            /// Create a stretching glyph (disc with hole), its center being the origin and extending in the x,y-plane with radius 1
            /// </summary>
            /// <param name="circle_resolution">Resolution of the circle (angle)</param>
            /// <param name="bending">Disc bending for pseudo-3D effect</param>
            /// <param name="hole_radius">Radius of the inner hole</param>
            /// <param name="offset">Offset in z-direction</param>
            /// <param name="strip_width">Width of the strip indicating the major eigenvectors</param>
            /// <param name="reference_width">Width of the reference circle</param>
            /// <param name="z_offset">Offset of the strips to prevent z-fighting</param>
            /// <returns>Template stretching glyph [disc, minimum eigenvalue strip, maximum eigenvalue strip, reference circle]</returns>
            stretching_glyph_t create_stretching_glyph_template(std::size_t circle_resolution, float_t bending, float_t hole_radius,
                float_t offset, float_t strip_width, float_t reference_width, float_t z_offset) const;

            /// <summary>
            /// "Instantiate" stretching glyphs at each interface position
            /// </summary>
            /// <param name="glyph_template">Template to instantiate</param>
            /// <param name="average_cell_size">Average cell size of the data</param>
            /// <param name="size_scalar">Scalar to modify the size relative to the average cell size</param>
            /// <param name="scalar">Scalar used for scaling the stretching factor</param>
            /// <param name="sharpness">Sharpness for superellipsoid</param>
            /// <param name="show_strips">Show strips indicating the major eigenvectors</param>
            /// <param name="show_reference">Show reference circle</param>
            void instantiate_stretching_glpyh(const stretching_glyph_t& glyph_template, float_t average_cell_size,
                float_t size_scalar, float_t scalar, float_t sharpness, bool show_strips, bool show_reference);

            /// <summary>
            /// Create a bending glyph (disc), its center being the origin and extending in the x,y-plane with radius 1
            /// </summary>
            /// <param name="circle_resolution">Resolution of the circle (angle)</param>
            /// <param name="polynomial_resolution">Resolution along the radius</param>
            /// <param name="offset">Offset in z-direction</param>
            /// <param name="strip_width">Width of the strip indicating the major eigenvectors</param>
            /// <param name="z_offset">Offset of the strips and reference circle to prevent z-fighting</param>
            /// <returns>Template bending glyph [disc, minimum eigenvalue strip, maximum eigenvalue strip]</returns>
            bending_glyph_t create_bending_glyph_template(std::size_t circle_resolution,
                std::size_t polynomial_resolution, float_t offset, float_t strip_width, float_t z_offset) const;

            /// <summary>
            /// "Instantiate" bending glyphs at each interface position
            /// </summary>
            /// <param name="glyph_template">Template to instantiate</param>
            /// <param name="average_cell_size">Average cell size of the data</param>
            /// <param name="size_scalar">Scalar to modify the size relative to the average cell size</param>
            /// <param name="scalar">Scalar to modify the bending</param>
            /// <param name="show_strips">Show strips indicating the major eigenvectors</param>
            void instantiate_bending_glyph(const bending_glyph_t& glyph_template,
                float_t average_cell_size, float_t size_scalar, float_t scalar, bool show_strips);

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

            data::polydata<float_t>* stretching_glyph_rings;
            data::polydata<float_t>* stretching_glyph_references;
            data::polydata<float_t>* stretching_glyph_strips;

            data::polydata<float_t>* bending_glyph_discs;
            data::polydata<float_t>* bending_glyph_strips;

            /// Create glyphs?
            bool velocity_glyph;
            bool stretching_glyph;
            bool bending_glyph;

            /// Time step for scaling the velocities
            float_t timestep;

            /// Properties of the velocity glyphs
            interface_deformation_glyph_aux::velocity_params_t<float_t> velocity_parameters;

            /// Properties of the stretching glyph
            interface_deformation_glyph_aux::stretching_params_t<float_t> stretching_parameters;

            /// Properties of the bending glyph
            interface_deformation_glyph_aux::bending_params_t<float_t> bending_parameters;
        };
    }
}

#include "tpf_module_interface_deformation_glyph.inl"
