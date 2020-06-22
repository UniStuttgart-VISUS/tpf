#include "tpf_module_stl_reader.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_geometric_object.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_math_defines.h"
#include "tpf/math/tpf_vector.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace tpf
{
    namespace modules
    {
       template <typename float_t, typename kernel_t>
        inline stl_reader<float_t, kernel_t>::stl_reader() { }

        template <typename float_t, typename kernel_t>
        inline std::string stl_reader<float_t, kernel_t>::get_name() const
        {
            return std::string("STL Reader");
        }

        template <typename float_t, typename kernel_t>
        inline void stl_reader<float_t, kernel_t>::set_algorithm_output(stl_reader_aux::polydata_or_buffer<float_t> output)
        {
            if (output.tag == output_t::POLYDATA)
            {
                this->triangles = output.polydata;
            }
            else
            {
                this->vertex_buffer = output.buffer;
            }
        }

        template <typename float_t, typename kernel_t>
        inline void stl_reader<float_t, kernel_t>::set_algorithm_parameters(const std::string& file_name)
        {
            this->file_name = file_name;
        }

        template <typename float_t, typename kernel_t>
        inline void stl_reader<float_t, kernel_t>::set_run_parameters(const std::size_t timestep, const data::extent_t& extent)
        { }

        template <typename float_t, typename kernel_t>
        inline void stl_reader<float_t, kernel_t>::set_additional_run_parameter(const bool calculate_normal)
        {
            this->calculate_normal = calculate_normal;
        }

        template <typename float_t, typename kernel_t>
        inline data_information stl_reader<float_t, kernel_t>::provide_information()
        {
            // Read STL header
            data_information info;

            std::ifstream ifs(this->file_name, std::ofstream::in | std::ofstream::binary);

            if (ifs.good())
            {
                uint32_t num_triangles;

                ifs.ignore(80 * sizeof(uint8_t));
                ifs.read(reinterpret_cast<char*>(&num_triangles), sizeof(uint32_t));

                info.datasets.push_back(std::make_pair("Surface Normals", static_cast<std::size_t>(num_triangles)));

                // Sanity check for file size
                ifs.ignore(std::numeric_limits<std::streamsize>::max());
                std::streamsize file_size = ifs.gcount();

                if (file_size != num_triangles * 50)
                {
                    throw std::runtime_error(__tpf_error_message("File size does not match the number of triangles."));
                }

                log::info_message(__tpf_info_message("Number of triangles from STL file: ", num_triangles));

                // Read file content for extent information
                ifs.seekg(80 * sizeof(uint8_t) + sizeof(uint32_t));

                math::vec3_t<float> min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
                math::vec3_t<float> max(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());

                for (std::size_t triangle_index = 0; triangle_index < static_cast<std::size_t>(num_triangles); ++triangle_index)
                {
                    math::vec3_t<float> vertex_1, vertex_2, vertex_3;

                    ifs.ignore(3 * sizeof(float));
                    ifs.read(reinterpret_cast<char*>(vertex_1.data()), 3 * sizeof(float));
                    ifs.read(reinterpret_cast<char*>(vertex_2.data()), 3 * sizeof(float));
                    ifs.read(reinterpret_cast<char*>(vertex_3.data()), 3 * sizeof(float));
                    ifs.ignore(sizeof(uint16_t));

                    min = minmax<true>(min, vertex_1, vertex_2, vertex_3);
                    max = minmax<false>(max, vertex_1, vertex_2, vertex_3);
                }

                info.bbox.resize(3);
                info.bbox[0].first = min[0];
                info.bbox[1].first = min[1];
                info.bbox[2].first = min[2];
                info.bbox[0].second = max[0];
                info.bbox[1].second = max[1];
                info.bbox[2].second = max[2];

                log::info_message(__tpf_info_message("Extent: [", min[0], ", ", min[1], ", ", min[2], "] x [", max[0], ", ", max[1], ", ", max[2], "]"));

                ifs.close();
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Error reading file header."));
            }

            return info;
        }

        template <typename float_t, typename kernel_t>
        inline void stl_reader<float_t, kernel_t>::run_algorithm()
        {
            if (this->triangles == nullptr && this->vertex_buffer == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Input must not be null."));
            }

            // Read STL file
            std::ifstream ifs(this->file_name, std::ofstream::in | std::ofstream::binary);

            if (ifs.good())
            {
                uint32_t num_triangles;

                ifs.ignore(80 * sizeof(uint8_t));
                ifs.read(reinterpret_cast<char*>(&num_triangles), sizeof(uint32_t));

#ifdef __tpf_sanity_checks
                bool normals_match = true;

#ifdef __tpf_detailed
                auto validity_ref = std::make_shared<data::array<unsigned char, 1>>(module_name::stl_reader::validity, static_cast<std::size_t>(num_triangles));
                auto& validity = *validity_ref;
#endif
#endif

                if (this->triangles != nullptr)
                {
                    std::vector<std::shared_ptr<geometry::geometric_object<float_t>>> triangles;
                    triangles.reserve(static_cast<std::size_t>(num_triangles));

                    auto normals_ref = std::make_shared<data::array<float_t, 3>>("Surface Normals", static_cast<std::size_t>(num_triangles));
                    auto& normals = *normals_ref;

                    for (std::size_t triangle_index = 0; triangle_index < static_cast<std::size_t>(num_triangles); ++triangle_index)
                    {
                        math::vec3_t<float> normal, vertex_1, vertex_2, vertex_3;

                        ifs.read(reinterpret_cast<char*>(normal.data()), 3 * sizeof(float));
                        ifs.read(reinterpret_cast<char*>(vertex_1.data()), 3 * sizeof(float));
                        ifs.read(reinterpret_cast<char*>(vertex_2.data()), 3 * sizeof(float));
                        ifs.read(reinterpret_cast<char*>(vertex_3.data()), 3 * sizeof(float));
                        ifs.ignore(sizeof(uint16_t));

#ifdef __tpf_sanity_checks
#ifndef __tpf_detailed
                        if (normals_match)
#endif
                        {
                            const math::vec3_t<float_t> calculated_normal = (vertex_2 - vertex_1).cross(vertex_3 - vertex_1).normalized();

                            const bool matches = calculated_normal.isApprox(normal) || math::calculate_angle(calculated_normal, normal) < (math::pi<float_t> / static_cast<float_t>(2.0));

                            normals_match &= matches;

#ifdef __tpf_detailed
                            validity(triangle_index) = static_cast<unsigned char>(matches ? 1 : 0);
#endif
                        }
#endif

                        if (this->calculate_normal)
                        {
                            normals(triangle_index) = (vertex_2 - vertex_1).cross(vertex_3 - vertex_1).normalized();
                        }
                        else
                        {
                            normals(triangle_index) = copy_vector<float_t>(normal);
                        }

                        geometry::point<float_t, kernel_t> point_1(copy_vector<float_t>(vertex_1));
                        geometry::point<float_t, kernel_t> point_2(copy_vector<float_t>(vertex_2));
                        geometry::point<float_t, kernel_t> point_3(copy_vector<float_t>(vertex_3));

                        triangles.push_back(std::make_shared<geometry::triangle<float_t, kernel_t>>(point_1, point_2, point_3));
                    }

                    this->triangles->insert(triangles);
                    this->triangles->add(normals_ref, data::topology_t::CELL_DATA);
                }
                else
                {
                    auto& vertex_buffer = *this->vertex_buffer;

                    const std::size_t header_block_size = sizeof(uint32_t);
                    const std::size_t data_block_size = 9 * num_triangles * sizeof(float);

                    const std::size_t data_offset_1 = 2 * header_block_size;
                    const std::size_t data_offset_2 = 2 * header_block_size + data_block_size;

                    vertex_buffer.resize(2 * data_block_size + 2 * header_block_size);

                    reinterpret_cast<uint32_t&>(vertex_buffer[0 * header_block_size]) = static_cast<uint32_t>(data_offset_1);
                    reinterpret_cast<uint32_t&>(vertex_buffer[1 * header_block_size]) = static_cast<uint32_t>(data_offset_2);

                    std::size_t vertex_position = data_offset_1;
                    std::size_t normal_position = data_offset_2;

                    for (std::size_t triangle_index = 0; triangle_index < static_cast<std::size_t>(num_triangles); ++triangle_index)
                    {
                        ifs.read(reinterpret_cast<char*>(&vertex_buffer[normal_position]), 3 * sizeof(float));
                        ifs.read(reinterpret_cast<char*>(&vertex_buffer[vertex_position]), 9 * sizeof(float));
                        ifs.ignore(sizeof(uint16_t));

                        std::memcpy(&vertex_buffer[normal_position + 3 * sizeof(float)], &vertex_buffer[normal_position], 3 * sizeof(float));
                        std::memcpy(&vertex_buffer[normal_position + 6 * sizeof(float)], &vertex_buffer[normal_position], 3 * sizeof(float));

#ifdef __tpf_sanity_checks
                        if (normals_match || this->calculate_normal)
#else
                        if (this->calculate_normal)
#endif
                        {

#ifdef __tpf_sanity_checks
                            const math::vec3_t<float> normal(
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 0 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 1 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 2 * sizeof(float)]));
#endif

                            const math::vec3_t<float> vertex_1(
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 0 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 1 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 2 * sizeof(float)]));

                            const math::vec3_t<float> vertex_2(
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 3 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 4 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 5 * sizeof(float)]));

                            const math::vec3_t<float> vertex_3(
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 6 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 7 * sizeof(float)]),
                                reinterpret_cast<float&>(vertex_buffer[vertex_position + 8 * sizeof(float)]));

                            const math::vec3_t<float> calculated_normal = (vertex_2 - vertex_1).cross(vertex_3 - vertex_1).normalized();
#ifdef __tpf_sanity_checks
                            normals_match &= calculated_normal.isApprox(normal) || math::calculate_angle(calculated_normal, normal) < (math::pi<float> / 2.0f);
#endif

                            if (this->calculate_normal)
                            {
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 0 * sizeof(float)]) = calculated_normal[0];
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 1 * sizeof(float)]) = calculated_normal[1];
                                reinterpret_cast<float&>(vertex_buffer[normal_position + 2 * sizeof(float)]) = calculated_normal[2];
                            }
                        }

                        vertex_position += 9 * sizeof(float);
                        normal_position += 9 * sizeof(float);
                    }
                }

#ifdef __tpf_sanity_checks
                if (!normals_match)
                {
                    log::warning_message(__tpf_warning_message("Normals do not match the vertex orientation"));
                }

#ifdef __tpf_detailed
                if (this->triangles != nullptr)
                {
                    this->triangles->add(validity_ref, data::topology_t::CELL_DATA);
                }
#endif
#endif

                ifs.close();
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Error reading data."));
            }
        }

        template <typename float_t, typename kernel_t>
        template <bool min, typename... T>
        inline math::vec3_t<float> stl_reader<float_t, kernel_t>::minmax(const math::vec3_t<float>& first, const math::vec3_t<float>& second, const T&... rest) const
        {
            math::vec3_t<float_t> min_or_max;

            min_or_max = minmax<min>(first, second);
            min_or_max = minmax<min>(min_or_max, rest...);

            return min_or_max;
        }

        template <typename float_t, typename kernel_t>
        template <bool min>
        inline math::vec3_t<float> stl_reader<float_t, kernel_t>::minmax(const math::vec3_t<float>& first, const math::vec3_t<float>& second) const
        {
            math::vec3_t<float_t> min_or_max;

            if (min)
            {
                min_or_max[0] = std::min(first[0], second[0]);
                min_or_max[1] = std::min(first[1], second[1]);
                min_or_max[2] = std::min(first[2], second[2]);
            }
            else
            {
                min_or_max[0] = std::max(first[0], second[0]);
                min_or_max[1] = std::max(first[1], second[1]);
                min_or_max[2] = std::max(first[2], second[2]);
            }

            return min_or_max;
        }
    }
}
