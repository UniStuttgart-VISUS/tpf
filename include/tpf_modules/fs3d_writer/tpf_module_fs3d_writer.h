#pragma once

#include "tpf_fs3d_data_writer.h"

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_parameters.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"

#include <memory>
#include <string>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module for reading FS3D list files
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <template name="Components">Number of components of the data to write</template>
        template <typename float_t, int Components>
        class fs3d_writer : public module_base<
            interface_input<const data::grid<float_t, float_t, 3, Components>&>,
            interface_parameters<const std::string&>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, Components>&>;
            using parameters_t = interface_parameters<const std::string&>;

            using base_t = module_base<input_t, parameters_t>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            fs3d_writer();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="grid">Input grid</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, Components>& grid) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="file_name">File name</param>
            virtual void set_algorithm_parameters(const std::string& file_name) override;

        public:
            /// <summary>
            /// Set parameters for module run
            /// </summary>
            /// <param name="timestep">Time step</param>
            /// <param name="global_extent">Global extent</param>
            void set_run_parameters(std::size_t timestep, const data::extent_t& global_extent, double time = 0.0,
                std::string name = "Data", std::string unit = "[-]", std::string grid_unit = "[cm]");

        protected:
            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// Fractions
            const data::grid<float_t, float_t, 3, Components>* grid;

            /// File names
            std::string file_name;
            std::string grid_file_name;

            /// Time step
            std::size_t timestep;
            double time;

            /// Extent
            data::extent_t global_extent;

            /// Measure and naming information
            std::string name, unit, grid_unit;

            /// Data reader
            std::unique_ptr<fs3d_writer_aux::fs3d_data_writer<float_t, float_t>> writer;

            /// Status of run parameter
            bool run_parameter_set;
        };
    }
}

#include "tpf_module_fs3d_writer.inl"
