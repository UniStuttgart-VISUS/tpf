#pragma once

#include "tpf_module_base.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include <string>
#include <utility>
#include <vector>

#include "tpf/data/tpf_data_information.h"

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Information on spatial and temporal extents
        /// </summary>
        struct data_information
        {
            /// Extent
            data::extent_t extent;

            /// Bounding box
            data::bbox_t<double> bbox;

            /// Time steps
            std::vector<double> timesteps;

            /// Data sets (name and number of components or objects)
            std::vector<std::pair<std::string, std::size_t>> datasets;
        };

        /// <summary>
        /// Reader module base class
        /// </summary>
        template <typename Output>
        class reader_base : public module_base<
            interface_output<Output>,
            interface_parameters<const std::string&>>
        {
            using output_t = interface_output<Output>;
            using parameters_t = interface_parameters<const std::string&>;

            using base_t = module_base<output_t, parameters_t>;

        public:
            /// <summary>
            /// Constructor
            /// </summary>
            reader_base();

            /// <summary>
            /// Set parameter for running the module
            /// </summary>
            /// <param name="timestep">Time step</param>
            /// <param name="extent">Extent</param>
            virtual void set_run_parameters(std::size_t timestep = 0, const data::extent_t& extent = data::extent_t(0)) = 0;

        protected:
            /// <summary>
            /// Provide information prior to module run
            /// </summary>
            virtual data_information provide_information() = 0;
        };
    }
}

#include "tpf_module_reader_base.inl"
