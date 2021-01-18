#pragma once

#include "tpf_particle_seed.h"
#include "tpf_stream_trace.h"

#include "tpf/data/tpf_polydata.h"

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_callbacks.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/policies/tpf_interpolatable.h"

#include <string>
#include <tuple>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            /// <summary>
            /// Call back class for requesting the next time frame
            /// </summary>
            template <typename float_t, typename point_t>
            class request_frame_call_back
            {
            public:
                /// <summary>
                /// Returns the data for the next time step if possible
                /// </summary>
                /// <returns>[Time step delta, velocity, global velocity part, translation, rotation, validity]</returns>
                virtual std::tuple<float_t,
                    policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>*,
                    policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>*,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<bool(const point_t&)>> operator()() = 0;

                /// <summary>
                /// Reset input algorithms to their beginning state
                /// </summary>
                virtual void reset() = 0;
            };

            /// <summary>
            /// method_t of vector flow visualization
            /// </summary>
            enum class method_t
            {
                INVALID = -1,
                STREAM, STREAK, PATH
            };
        }

        /// <summary>
        /// Module for calculation of flow fields
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <template name="point_t">Point type</template>
        template <typename float_t, typename point_t>
        class flow_field : public module_base<
            interface_callbacks<flow_field_aux::request_frame_call_back<float_t, point_t>*>,
            interface_input<const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>&,
                const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>&,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                std::function<bool(const point_t&)>, const data::polydata<float_t>&>,
            interface_output<data::polydata<float_t>&>,
            interface_parameters<flow_field_aux::method_t, std::size_t, float_t>>
        {
        public:
            using callbacks_t = interface_callbacks<flow_field_aux::request_frame_call_back<float_t, point_t>*>;
            using input_t = interface_input<const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>&,
                const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>&,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                std::function<bool(const point_t&)>, const data::polydata<float_t>&>;
            using output_t = interface_output<data::polydata<float_t>&>;
            using parameters_t = interface_parameters<flow_field_aux::method_t, std::size_t, float_t>;

            using base_t = module_base<callbacks_t, input_t, output_t, parameters_t>;

            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            flow_field();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</param>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set callbacks
            /// </summary>
            /// <param name="next_time_frame_callback">Callback for requesting next time frame</param>
            virtual void set_algorithm_callbacks(flow_field_aux::request_frame_call_back<float_t, point_t>* next_time_frame_callback) override;

            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="velocities">Input velocities</param>
            /// <param name="global_velocity_parts">Input global velocity parts</param
            /// <param name="get_translation">Function to get translation at particle position</param>
            /// <param name="get_rotation_axis">Function to get rotation at particle position</param>
            /// <param name="is_particle_valid">Function to check the validity of particles</param>
            /// <param name="seed">Input seed</param>
            virtual void set_algorithm_input(const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>& velocities,
                const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>& global_velocity_parts,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_translation,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_rotation_axis,
                std::function<bool(const point_t&)> is_particle_valid, const data::polydata<float_t>& seed) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="lines">Output flow lines ("ID (Advection)", "ID (Distribution)")</param>
            virtual void set_algorithm_output(data::polydata<float_t>& lines) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="method">method_t to use: 0 - streamlines, 1 - streaklines, 2 - pathlines</param>
            /// <param name="num_advections">Number of advections</param>
            /// <param name="timestep">Time step</param>
            virtual void set_algorithm_parameters(flow_field_aux::method_t method, std::size_t num_advections, float_t timestep) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Compute streamlines
            /// </summary>
            /// <param name="seed">Particle seed</param>
            void compute_streamlines(flow_field_aux::particle_seed<float_t>&& seed);

            /// <summary>
            /// Compute pathlines
            /// </summary>
            /// <param name="seed">Particle seed</param>
            void compute_pathlines(flow_field_aux::particle_seed<float_t>&& seed);

            /// <summary>
            /// Compute streaklines
            /// </summary>
            /// <param name="seed">Particle seed</param>
            void compute_streaklines(flow_field_aux::particle_seed<float_t>&& seed);

            /// <summary>
            /// Create and save line segments for all particle advections
            /// </summary>
            /// <param name="particles">Particle traces</param>
            /// <param name="num_advections">Number of the current advection step</param>
            /// <param name="inverse">Invert direction of lines</param>
            void create_lines(const flow_field_aux::stream_trace<float_t>& particles, std::size_t num_advections, bool inverse = false);

            /// Velocities
            const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>* velocities;
            const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>* global_velocity_parts;

            /// Callback functions
            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_translation;
            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_rotation_axis;
            std::function<bool(const point_t&)> is_particle_valid;

            /// Seed
            const data::polydata<float_t>* seed;

            /// Integration lines
            data::polydata<float_t>* lines;

            /// method_t to use
            flow_field_aux::method_t method;

            /// Number of advections
            std::size_t num_advections;

            /// Time step
            float_t timestep;

            /// Callback for requesting next time frame
            flow_field_aux::request_frame_call_back<float_t, point_t>* next_time_frame_callback;
        };
    }
}

#include "tpf_module_flow_field.inl"
