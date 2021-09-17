#pragma once

#include "tpf_particle_seed.h"
#include "tpf_stream_trace.h"

#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_callbacks.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/policies/tpf_interpolatable.h"

#include <functional>
#include <string>
#include <tuple>
#include <vector>

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
                /// <returns>[Time step delta, sample time, velocity, translation, angular velocity, barycenter,
                ///  initial translation, initial angular velocity, validity,
                ///  fields to interpolate and store at the particle positions]</returns>
                virtual std::tuple<float_t, float_t,
                    policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>*,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)>,
                    std::function<bool(const point_t&)>,
                    std::vector<std::tuple<std::string, std::size_t, policies::interpolatable_base<point_t>*>>> operator()() = 0;

                /// <summary>
                /// Reset input algorithms to their beginning state
                /// </summary>
                virtual void reset() = 0;
            };

            /// <summary>
            /// Method of flow visualization
            /// </summary>
            enum class method_t
            {
                INVALID = -1,
                STREAM, STREAK, PATH
            };

            /// <summary>
            /// Dynamic vs. static frame of reference
            /// </summary>
            enum class time_dependency_t
            {
                dynamic,
                static_
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
            interface_input<const data::polydata<float_t>&>,
            interface_output<data::polydata<float_t>&>,
            interface_parameters<flow_field_aux::method_t, std::size_t,
                flow_field_aux::time_dependency_t, bool, bool>>
        {
        public:
            using callbacks_t = interface_callbacks<flow_field_aux::request_frame_call_back<float_t, point_t>*>;
            using input_t = interface_input<const data::polydata<float_t>&>;
            using output_t = interface_output<data::polydata<float_t>&>;
            using parameters_t = interface_parameters<flow_field_aux::method_t, std::size_t,
                flow_field_aux::time_dependency_t, bool, bool>;

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
            /// <param name="seed">Input seed</param>
            virtual void set_algorithm_input(const data::polydata<float_t>& seed) override;

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
            /// <param name="time_dependency">Dynamic vs. static frame of reference</param>
            /// <param name="keep_translation">Keep translational velocity part</param>
            /// <param name="keep_rotation">Keep rotational velocity part</param>
            virtual void set_algorithm_parameters(flow_field_aux::method_t method, std::size_t num_advections,
                flow_field_aux::time_dependency_t time_dependency, bool keep_translation, bool keep_rotation) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Get global velocity parts for the particle
            /// </summary>
            /// <param name="particle">Particle position</param>
            /// <param name="get_translation">Function to get translational velocity part</param>
            /// <param name="get_angular_velocity">Function to get angular velocity</param>
            /// <param name="get_barycenter">Function to get droplet barycenter</param>
            /// <param name="get_translation">Function to get initial translational velocity part</param>
            /// <param name="get_angular_velocity">Function to get initial angular velocity</param>
            /// <returns>Global velocity part</returns>
            Eigen::Matrix<float_t, 3, 1> get_global_velocity(const point_t& particle,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_translation,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_angular_velocity,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_barycenter,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_initial_translation,
                std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_initial_angular_velocity) const;

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
            /// <param name="time">Time information</param>
            /// <param name="inverse">Invert direction of lines</param>
            void create_lines(const flow_field_aux::stream_trace<float_t>& particles, std::size_t num_advections,
                const std::vector<double>& time, bool inverse = false);

            /// Seed
            const data::polydata<float_t>* seed;

            /// Integration lines
            data::polydata<float_t>* lines;

            /// Method to use
            flow_field_aux::method_t method;

            /// Number of advections
            std::size_t num_advections;

            /// Dynamic vs. static frame of reference
            flow_field_aux::time_dependency_t time_dependency;

            /// Keep velocity parts?
            bool keep_translation;
            bool keep_rotation;

            /// Callback for requesting next time frame
            flow_field_aux::request_frame_call_back<float_t, point_t>* next_time_frame_callback;
        };
    }
}

#include "tpf_module_flow_field.inl"
