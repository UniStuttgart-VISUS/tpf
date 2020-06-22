#pragma once

#include "tpf/module/tpf_module_reader_base.h"

#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_geometric_object.h"

#include "tpf/math/tpf_vector.h"

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace stl_reader_aux
        {
            template <typename float_t>
            class polydata_or_buffer
            {
            public:
                enum
                {
                    POLYDATA, BUFFER
                } tag;

                union
                {
                    data::polydata<float_t>* polydata;
                    std::vector<uint8_t>* buffer;
                };
            };
        }

        /// <summary>
        /// Module for reading FS3D list files
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel for triangle objects</template>
        template <typename float_t, typename kernel_t = geometry::default_kernel_t>
        class stl_reader : public reader_base<stl_reader_aux::polydata_or_buffer<float_t>>
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            stl_reader();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            virtual std::string get_name() const override;

        protected:
            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="output">Output triangles or vertex buffer</param>
            virtual void set_algorithm_output(stl_reader_aux::polydata_or_buffer<float_t> output) override;

            /// <summary>
            /// Set parameters
            /// </summary>
            /// <param name="file_name">File name</param>
            virtual void set_algorithm_parameters(const std::string& file_name) override;

            /// <summary>
            /// Set parameters for module run
            /// </summary>
            /// <param name="timestep">Time step</param>
            /// <param name="extent">Extent</param>
            virtual void set_run_parameters(std::size_t timestep = 0, const data::extent_t& extent = data::extent_t(0)) override;

        public:
            /// <summary>
            /// Set additional parameters for module run
            /// </summary>
            /// <param name="calculate_normal">Calculate normal instead of reading it from file</param>
            void set_additional_run_parameter(bool calculate_normal);

            /// <summary>
            /// Provide information prior to module run
            /// </summary>
            data_information provide_information() override;

        protected:
            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Return the provided vector unmodified
            /// </summary>
            /// <template name="OT">Same data type as NT</template>
            /// <param name="input">Input vector</param>
            /// <returns>Original input vector</returns>
            template <typename OT>
            math::vec3_t<float_t>& copy_vector(math::vec3_t<OT>& input, typename std::enable_if<std::is_same<OT, float_t>::value>::type* = nullptr) const
            {
                return input;
            }

            /// <summary>
            /// Return the provided vector unmodified
            /// </summary>
            /// <template name="OT">Same data type as NT</template>
            /// <param name="input">Input vector</param>
            /// <returns>Original input vector</returns>
            template <typename OT>
            const math::vec3_t<float_t>& copy_vector(const math::vec3_t<OT>& input, typename std::enable_if<std::is_same<OT, float_t>::value>::type* = nullptr) const
            {
                return input;
            }

            /// <summary>
            /// Copy the vector to a new type
            /// </summary>
            /// <template name="OT">Original data type</template>
            /// <param name="input">Input vector</param>
            /// <returns>Copied vector</returns>
            template <typename OT>
            math::vec3_t<float_t> copy_vector(const math::vec3_t<OT>& input, typename std::enable_if<!std::is_same<OT, float_t>::value>::type* = nullptr) const
            {
                return math::vec3_t<float_t>(static_cast<float_t>(input[0]), static_cast<float_t>(input[1]), static_cast<float_t>(input[2]));
            }

            /// <summary>
            /// Calculate the minimum or maximum from a vector
            /// </summary>
            /// <param name="first">First vector</param>
            /// <param name="second">Second vector</param>
            /// <param name="...rest">More vectors</param>
            /// <returns>Minimum or maximum</returns>
            template <bool min, typename... T>
            math::vec3_t<float> minmax(const math::vec3_t<float>& first, const math::vec3_t<float>& second, const T&... rest) const;

            template <bool min>
            math::vec3_t<float> minmax(const math::vec3_t<float>& first, const math::vec3_t<float>& second) const;

            /// Triangles
            data::polydata<float_t>* triangles;

            /// Vertex buffer
            std::vector<uint8_t>* vertex_buffer;

            /// File names
            std::string file_name;

            /// Calculate normal or use the one from file
            bool calculate_normal;
        };
    }
}

#include "tpf_module_stl_reader.inl"
