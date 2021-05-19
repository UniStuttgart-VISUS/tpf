#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_callbacks.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_polydata.h"

#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace dynamic_droplets_aux
        {
            /// <summary>
            /// Callback class for requesting the next time frame
            /// </summary>
            /// <template name="float_t">Floating point type</template>
            template <typename float_t>
            class request_frame_call_back
            {
            public:
                /// <summary>
                /// Returns the data for the next time step if possible
                /// </summary>
                /// <returns>[Time step delta, droplets]</returns>
                virtual std::pair<float_t, data::polydata<float_t>> operator()() = 0;

                /// <summary>
                /// Reset input algorithms to their beginning state
                /// </summary>
                virtual void reset() = 0;
            };
        }

        /// <summary>
        /// Module for calculation of droplet context
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class dynamic_droplets : public module_base<
            interface_callbacks<dynamic_droplets_aux::request_frame_call_back<float_t>*>,
            interface_input<const data::polydata<float_t>&>,
            interface_output<
                opt_arg<data::polydata<float_t>>,
                opt_arg<data::polydata<float_t>>,
                opt_arg<data::polydata<float_t>>>,
            interface_parameters<std::size_t, float_t, float_t, bool, std::string, std::string, std::string>>
        {
        public:
            using callbacks_t = interface_callbacks<dynamic_droplets_aux::request_frame_call_back<float_t>*>;
            using input_t = interface_input<const data::polydata<float_t>&>;
            using output_t = interface_output<
                opt_arg<data::polydata<float_t>>,
                opt_arg<data::polydata<float_t>>,
                opt_arg<data::polydata<float_t>>>;
            using parameters_t = interface_parameters<std::size_t, float_t, float_t,
                bool, std::string, std::string, std::string>;

            using base_t = module_base<callbacks_t, input_t, output_t, parameters_t>;

            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            dynamic_droplets();

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
            virtual void set_algorithm_callbacks(dynamic_droplets_aux::request_frame_call_back<float_t>* next_time_frame_callback) override;

            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="droplets">Input droplets</param>
            virtual void set_algorithm_input(const tpf::data::polydata<float_t>& droplets) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="paths">Output translation paths</param>
            /// <param name="axes">Output rotation axes</param>
            /// <param name="ribbons">Output rotation ribbons</param>
            virtual void set_algorithm_output(opt_arg<data::polydata<float_t>> paths,
                opt_arg<data::polydata<float_t>> axes, opt_arg<data::polydata<float_t>> ribbons) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="num_timesteps">Number of time steps</param>
            /// <param name="timestep">Time step</param>
            /// <param name="axis_scale">Scale factor for rotation axes</param>
            /// <param name="axis_translation">Translate axis back to its origin</param>
            virtual void set_algorithm_parameters(std::size_t num_timesteps, float_t timestep, float_t axis_scale,
                bool axis_translation, std::string translation_name, std::string rotation_name, std::string radius_name) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// Droplet information
            struct droplet_t
            {
                Eigen::Matrix<float_t, 3, 1> position;
                Eigen::Matrix<float_t, 3, 1> translation;
                Eigen::Matrix<float_t, 3, 1> rotation;
                float_t radius;
            };

            /// <summary>
            /// Create and save translation paths
            /// </summary>
            /// <param name="droplets">Droplets storing the translation</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_paths(const std::vector<tpf::data::polydata<float_t>>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save rotation axes
            /// </summary>
            /// <param name="droplets">Droplets storing the rotation axes</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_axes(const std::vector<tpf::data::polydata<float_t>>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save rotation ribbons
            /// </summary>
            /// <param name="droplets">Droplets storing the rotation axes</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_ribbons(const std::vector<tpf::data::polydata<float_t>>& droplets, const std::vector<float_t>& timesteps);

            /// Droplets
            const tpf::data::polydata<float_t>* droplets;

            /// Translation paths
            data::polydata<float_t>* paths;

            /// Rotation axes and ribbons
            data::polydata<float_t>* axes;
            data::polydata<float_t>* ribbons;

            /// Number of time steps
            std::size_t num_timesteps;

            /// Time step
            float_t timestep;

            /// Scale factor for rotation axes
            float_t axis_scale;

            /// Translate axes back to their origin
            bool axis_translation;

            /// Array names
            std::string translation_name;
            std::string rotation_name;
            std::string radius_name;

            /// Callback for requesting next time frame
            dynamic_droplets_aux::request_frame_call_back<float_t>* next_time_frame_callback;
        };
    }
}

#include "tpf_module_dynamic_droplets.inl"
