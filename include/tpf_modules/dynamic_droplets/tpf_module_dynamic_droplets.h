#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_callbacks.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include <functional>
#include <optional>
#include <string>
#include <tuple>
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
                /// <returns>[Time step delta, droplets, droplet IDs]</returns>
                virtual std::tuple<float_t, data::polydata<float_t>, data::grid<long long, float_t, 3, 1>> operator()() = 0;

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
            interface_input<
                const data::polydata<float_t>&,
                const data::grid<long long, float_t, 3, 1>&>,
            interface_output<
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&>,
            interface_parameters<std::size_t, float_t, bool, float_t, bool, float_t, bool, std::string, std::string, std::string>>
        {
        public:
            using callbacks_t = interface_callbacks<dynamic_droplets_aux::request_frame_call_back<float_t>*>;
            using input_t = interface_input<
                const data::polydata<float_t>&,
                const data::grid<long long, float_t, 3, 1>&>;
            using output_t = interface_output<
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&,
                data::polydata<float_t>&>;
            using parameters_t = interface_parameters<std::size_t, float_t, bool, float_t, bool, float_t,
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
            /// <param name="droplet_ids">Input droplet IDs</param>
            virtual void set_algorithm_input(const tpf::data::polydata<float_t>& droplets,
                const data::grid<long long, float_t, 3, 1>& droplet_ids) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="tracks">Output droplet tracking information</param>
            /// <param name="summary">Output droplet summary</param>
            /// <param name="paths">Output translation paths</param>
            /// <param name="axes">Output rotation axes</param>
            /// <param name="ribbons">Output rotation ribbons</param>
            /// <param name="rotation_paths">Output rotation paths</param>
            /// <param name="coordinate_axes">Output transforming coordinate axes</param>
            virtual void set_algorithm_output(data::polydata<float_t>& tracks, data::polydata<float_t>& summary,
                data::polydata<float_t>& paths, data::polydata<float_t>& axes, data::polydata<float_t>& ribbons,
                data::polydata<float_t>& rotation_paths, data::polydata<float_t>& coordinate_axes) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="num_timesteps">Number of time steps</param>
            /// <param name="timestep">Time step</param>
            /// <param name="static_frame_of_reference">Static frame of reference instead of dynamic</param>
            /// <param name="ribbon_scale">Scale factor for ribbons</param>
            /// <param name="fix_axis_size">Fix axis length to circumsphere</param>
            /// <param name="axis_scale">Scale factor for rotation axes</param>
            /// <param name="axis_translation">Translate axis back to its origin</param>
            /// <param name="translation_name">Name of the array containing translation information</param>
            /// <param name="rotation_name">Name of the array containing rotation information</param>
            /// <param name="radius_name">Name of the array containing droplet radius information</param>
            virtual void set_algorithm_parameters(std::size_t num_timesteps, float_t timestep, bool static_frame_of_reference,
                float_t ribbon_scale, bool fix_axis_size, float_t axis_scale, bool axis_translation, std::string translation_name,
                std::string rotation_name, std::string radius_name) override;

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

            /// Tracking information
            enum class topology_t
            {
                none,
                breakup,
                collision
            };

            using droplet_trace_t = std::pair<std::vector<droplet_t>, topology_t>;

            /// <summary>
            /// Track droplets over time and return this trace
            /// </summary>
            /// <returns>Droplets over time</returns>
            std::pair<std::vector<droplet_trace_t>, std::vector<float_t>> track_droplets() const;

            /// <summary>
            /// Create static droplet summary
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_summary(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save translation paths
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_paths(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save rotation axes
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_axes(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save rotation ribbons
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_ribbons(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save rotation paths
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_rotation_paths(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// <summary>
            /// Create and save transforming coordinate axes
            /// </summary>
            /// <param name="droplets">Droplet information over time</param>
            /// <param name="timesteps">Time step sizes</param>
            void create_coordinate_axes(const std::vector<droplet_trace_t>& droplets, const std::vector<float_t>& timesteps);

            /// Droplets
            const tpf::data::polydata<float_t>* droplets;
            const tpf::data::grid<long long, float_t, 3, 1>* droplet_ids;

            /// Droplet tracks
            data::polydata<float_t>* tracks;

            /// Droplet summary
            data::polydata<float_t>* summary;

            /// Translation paths
            data::polydata<float_t>* paths;

            /// Rotation axes, ribbons and particles
            data::polydata<float_t>* axes;
            data::polydata<float_t>* ribbons;
            data::polydata<float_t>* rotation_paths;

            /// Coordinate system axes
            data::polydata<float_t>* coordinate_axes;

            /// Number of time steps
            std::size_t num_timesteps;

            /// Time step
            float_t timestep;

            /// Static frame of reference?
            bool static_frame_of_reference;

            /// Scale factor for ribbons
            float_t ribbon_scale;

            /// Scale factor for rotation axes
            bool fix_axis_size;
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
