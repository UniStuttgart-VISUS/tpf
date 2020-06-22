#include "tpf_fs3d_data_writer.h"

#include "../fs3d_reader/tpf_fs3d_data.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"

#include "tpf/performance/tpf_timer.h"

#include "tpf/exception/tpf_not_implemented_exception.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_writer_aux
        {
            template <typename value_t, typename point_t>
            inline fs3d_data_writer<value_t, point_t>::fs3d_data_writer(const std::string& data_filename, const std::string& grid_filename)
                : data_filename(data_filename), grid_filename(grid_filename)
            { }

            template <typename value_t, typename point_t>
            inline void fs3d_data_writer<value_t, point_t>::set_global_extent(const data::extent_t& global_extent)
            {
                this->global_extent = global_extent;
            }

            template <typename value_t, typename point_t>
            template <int Components>
            inline void fs3d_data_writer<value_t, point_t>::write_timestep(const std::size_t timestep, const double time, const std::string& name,
                const std::string& unit, const std::string& grid_unit, const data::grid<value_t, point_t, 3, Components>& grid)
            {
                const auto& mpi = mpi::get_instance();

                this->subextent = grid.get_extent();

                // Write data information
                fs3d_reader_aux::data_header data_header;
                data_header.components = Components;
                data_header.identifier = 0;
                data_header.name = name;
                data_header.time_id = static_cast<int>(timestep);
                data_header.time_step = time;
                data_header.unit = unit;
                data_header.x_resolution = static_cast<int>(grid.get_extent()[0].second - grid.get_extent()[0].first + 1);
                data_header.y_resolution = static_cast<int>(grid.get_extent()[1].second - grid.get_extent()[1].first + 1);
                data_header.z_resolution = static_cast<int>(grid.get_extent()[2].second - grid.get_extent()[2].first + 1);

                write_data_header(timestep, data_header);

                // Write data
                if (mpi.get_num_processes() == 1)
                {
                    write_data<Components>(timestep, grid.get_data());
                }
                else
                {
                    write_data_block<Components>(timestep, grid.get_data());
                }

                // If necessary, write grid file
                if (!this->grid_filename.empty())
                {
                    // Write grid information
                    fs3d_reader_aux::grid_header grid_header;
                    grid_header.unit = grid_unit;
                    grid_header.x_resolution = static_cast<int>(grid.get_extent()[0].second - grid.get_extent()[0].first + 1);
                    grid_header.y_resolution = static_cast<int>(grid.get_extent()[1].second - grid.get_extent()[1].first + 1);
                    grid_header.z_resolution = static_cast<int>(grid.get_extent()[2].second - grid.get_extent()[2].first + 1);

                    write_grid_header(grid_header);

                    // Write grid
                    if (mpi.get_num_processes() == 1)
                    {
                        write_grid(grid.get_grid_information());
                    }
                    else
                    {
                        throw tpf::exception::not_implemented_exception();
                    }
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_writer<value_t, point_t>::write_grid_header(const fs3d_reader_aux::grid_header& header) const
            {
                // Write output header
                const std::size_t max_string = 80;

                std::ofstream ofs(this->grid_filename, std::ofstream::out | std::ofstream::binary);

                if (ofs.good())
                {
                    std::vector<char> empty(max_string, '\0');

                    ofs.write(header.unit.c_str(), std::min(header.unit.size(), max_string));
                    ofs.write(empty.data(), max_string - std::min(header.unit.size(), max_string));

                    ofs.write(reinterpret_cast<const char*>(&header.x_resolution), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.y_resolution), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.z_resolution), sizeof(int));

                    ofs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error writing grid file header!"));
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_writer<value_t, point_t>::write_grid(const data::grid_information<point_t>& grid) const
            {
                // Write data
                std::ofstream ofs(this->grid_filename, std::ofstream::app | std::ofstream::binary);

                // Create output
                if (ofs.good())
                {
                    // Write for each dimension
                    if (std::is_same<typename fs3d_reader_aux::grid_header::value_type, point_t>::value)
                    {
                        // Write block
                        ofs.write((char*)grid.cell_coordinates[0].data(), grid.cell_coordinates[0].size() * sizeof(typename fs3d_reader_aux::grid_header::value_type));
                        ofs.write((char*)grid.cell_coordinates[1].data(), grid.cell_coordinates[1].size() * sizeof(typename fs3d_reader_aux::grid_header::value_type));
                        ofs.write((char*)grid.cell_coordinates[2].data(), grid.cell_coordinates[2].size() * sizeof(typename fs3d_reader_aux::grid_header::value_type));
                    }
                    else
                    {
                        for (std::size_t d = 0; d < 3; ++d)
                        {
                            // Fill buffer
                            std::vector<typename fs3d_reader_aux::grid_header::value_type> buffer(grid.cell_coordinates[d].size());

                            for (std::size_t i = 0; i < buffer.size(); ++i)
                            {
                                buffer[i] = static_cast<typename fs3d_reader_aux::grid_header::value_type>(grid.cell_coordinates[d][i]);
                            }

                            // Write buffer to file
                            ofs.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(typename fs3d_reader_aux::grid_header::value_type));
                        }
                    }
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error writing grid file!"));
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_writer<value_t, point_t>::write_data_header(const std::size_t timestep, const fs3d_reader_aux::data_header& header) const
            {
                // Create file name with time step suffix
                const std::string leading_zeroes("0000000000");
                std::string suffix("_");

                const std::string number = std::to_string(timestep);
                const std::string number_with_lead = leading_zeroes + number;
                suffix += number_with_lead.substr(std::min(number.length(), leading_zeroes.length()));

                const auto dot_pos = this->data_filename.find_last_of('.');
                const std::string file_name = this->data_filename.substr(0, dot_pos) + suffix + this->data_filename.substr(dot_pos);

                // Write output header
                const std::size_t max_string = 80;

                std::ofstream ofs(file_name, std::ofstream::out | std::ofstream::binary);

                if (ofs.good())
                {
                    std::vector<char> empty(max_string, '\0');

                    ofs.write(header.name.c_str(), std::min(header.name.size(), max_string));
                    ofs.write(empty.data(), max_string - std::min(header.name.size(), max_string));

                    ofs.write(header.unit.c_str(), std::min(header.unit.size(), max_string));
                    ofs.write(empty.data(), max_string - std::min(header.unit.size(), max_string));

                    ofs.write(reinterpret_cast<const char*>(&header.time_id), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.time_step), sizeof(double));
                    ofs.write(reinterpret_cast<const char*>(&header.x_resolution), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.y_resolution), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.z_resolution), sizeof(int));
                    ofs.write(reinterpret_cast<const char*>(&header.identifier), sizeof(int));

                    ofs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error writing data file header!"));
                }
            }

            template <typename value_t, typename point_t>
            template <int Components>
            inline void fs3d_data_writer<value_t, point_t>::write_data(const std::size_t timestep, const std::vector<value_t>& data) const
            {
                __tpf_timer_start(timer, "write_data", std::chrono::seconds);

                // Create file name with time step suffix
                const std::string leading_zeroes("0000000000");
                std::string suffix("_");

                const std::string number = std::to_string(timestep);
                const std::string number_with_lead = leading_zeroes + number;
                suffix += number_with_lead.substr(std::min(number.length(), leading_zeroes.length()));

                const auto dot_pos = this->data_filename.find_last_of('.');
                const std::string file_name = this->data_filename.substr(0, dot_pos) + suffix + this->data_filename.substr(dot_pos);

                // Open file stream
                std::ofstream ofs(file_name, std::ofstream::app | std::ofstream::binary);

                if (ofs.good())
                {
                    // Write for each component
                    if (std::is_same<typename fs3d_reader_aux::data_header::value_type, value_t>::value && Components == 1)
                    {
                        // Write block
                        ofs.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(typename fs3d_reader_aux::data_header::value_type));
                    }
                    else
                    {
                        std::vector<typename fs3d_reader_aux::data_header::value_type> buffer(data.size() / Components);

                        for (std::size_t c = 0; c < Components; ++c)
                        {
                            // Copy data to buffer
                            for (std::size_t i = 0; i < buffer.size(); ++i)
                            {
                                const std::size_t data_index = (i * Components) + c;

                                buffer[i] = static_cast<typename fs3d_reader_aux::data_header::value_type>(data[data_index]);
                            }

                            // Write block
                            ofs.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(typename fs3d_reader_aux::data_header::value_type));
                        }
                    }

                    ofs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error writing data file!"));
                }
            }

            template <typename value_t, typename point_t>
            template <int Components>
            inline void fs3d_data_writer<value_t, point_t>::write_data_block(const std::size_t timestep, const std::vector<value_t>& data) const
            {
                throw tpf::exception::not_implemented_exception();
            }
        }
    }
}
