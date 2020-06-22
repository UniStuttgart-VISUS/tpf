#pragma once

#include "tpf/module/tpf_module_base.h"

#include "tpf_fs3d_data.h"
#include "tpf_fs3d_data_reader.h"

#include <memory>
#include <string>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_reader_aux
        {
            template <typename float_t>
            class scalar_or_vector
            {
            public:
                enum
                {
                    SCALAR, VECTOR
                } tag;

                union
                {
                    data::grid<float_t, float_t, 3, 1>* scalar_field;
                    data::grid<float_t, float_t, 3, 3>* vector_field;
                };
            };
        }

        /// <summary>
        /// Module for reading FS3D list files
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class fs3d_reader : public reader_base<fs3d_reader_aux::scalar_or_vector<float_t>>
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            fs3d_reader();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="output">Output scalar field or vector field</param>
            virtual void set_algorithm_output(fs3d_reader_aux::scalar_or_vector<float_t> output) override;

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
            /// <param name="extent">Extent</param>
            virtual void set_run_parameters(std::size_t timestep, const data::extent_t& extent) override;

            /// <summary>
            /// Provide information prior to module run
            /// </summary>
            data_information provide_information() override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// Fractions
            data::grid<float_t, float_t, 3, 1>* fractions;

            /// Velocities
            data::grid<float_t, float_t, 3, 3>* velocities;

            /// File names
            std::string file_name;
            std::string grid_file_name;

            /// Time step
            std::size_t timestep;

            /// Extent
            data::extent_t extent;

            /// Data reader
            std::unique_ptr<fs3d_reader_aux::fs3d_data_reader<float_t, float_t>> reader;
        };
    }
}

#include "tpf_module_fs3d_reader.inl"
