#include "tpf_fs3d_data_reader.h"

#include "tpf_fs3d_data.h"

#include "tpf/module/tpf_module_reader_base.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"

#include "tpf/log/tpf_log.h"

#include "tpf/performance/tpf_timer.h"

#ifdef __tpf_use_mpiio
#include "tpf/utility/tpf_copy.h"
#include "tpf/utility/tpf_downcast.h"

#include "mpi.h"
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_reader_aux
        {
            template <typename value_t, typename point_t>
            inline fs3d_data_reader<value_t, point_t>::fs3d_data_reader(const std::string& list_filename, const std::string& grid_filename)
                : list_filename(list_filename), grid_filename(grid_filename), state(NO_INFORMATION)
            { }

            template <typename value_t, typename point_t>
            inline data_information fs3d_data_reader<value_t, point_t>::get_information()
            {
                // Read information
                read_information();

                // Return information
                data_information info;
                info.extent = this->extent;
                info.timesteps = this->timesteps;
                info.datasets.push_back(std::make_pair(this->dataset_name, this->num_components));

                return info;
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_reader<value_t, point_t>::set_subextent(const data::extent_t& extent)
            {
                if (this->state == NO_INFORMATION)
                {
                    throw std::runtime_error(__tpf_error_message("Global information must be set before setting subextents!"));
                }

                // Check for same sub extent and same number of ghost levels
                bool same_subextent = extent.size() == this->subextent.size();

                for (std::size_t d = 0; d < std::min(extent.size(), this->subextent.size()); ++d)
                {
                    same_subextent &= extent[d].first == this->subextent[d].first && extent[d].second == this->subextent[d].second;
                }

                // Only if there are changes, reload the grid
                if (!same_subextent)
                {
                    this->subextent = extent;

                    this->state = EXTENT_INFORMATION;
                }
            }

            template <typename value_t, typename point_t>
            template <std::size_t components>
            inline data::grid<value_t, point_t, 3, components> fs3d_data_reader<value_t, point_t>::read_timestep(const std::size_t timestep)
            {
                const auto& mpi = mpi::get_instance();

                // Read information
                read_information();

                // Check time step
                if (timestep >= this->data_filenames.size())
                {
                    throw std::runtime_error(__tpf_error_message("Time step is out of bounds!"));
                }

                // Check for same (sub) extent
                bool same_subextent = this->extent.size() == this->subextent.size();

                for (std::size_t d = 0; d < std::min(this->extent.size(), this->subextent.size()); ++d)
                {
                    same_subextent &= this->extent[d].first == this->subextent[d].first && this->extent[d].second == this->subextent[d].second;
                }

                // Read data information
                const auto data_header = read_data_header(timestep);

                // Sanity check
                if (data_header.x_resolution != this->header.x_resolution ||
                    data_header.y_resolution != this->header.y_resolution ||
                    data_header.z_resolution != this->header.z_resolution)
                {
                    throw std::runtime_error(__tpf_error_message("Grid resolutions do not match!"));
                }

                // Read data
                std::vector<value_t> data;

                if (mpi.get_num_processes() == 1 && same_subextent)
                {
                    data = read_data<components>(timestep);
                }
                else
                {
                    data = read_data_block<components>(timestep);
                }

                // Create grid
                return data::grid<value_t, point_t, 3, components>(data_header.name + data_header.unit, this->subextent, std::move(data),
                    this->cell_coordinates, this->node_coordinates, this->cell_sizes);
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_reader<value_t, point_t>::read_information()
            {
                if (this->state != ALL_INFORMATION)
                {
                    const auto& mpi = mpi::get_instance();

                    if (this->state == GLOBAL_INFORMATION)
                    {
                        throw std::runtime_error(__tpf_error_message("Subextent must be set before reading grid information!"));
                    }

                    if (this->state == NO_INFORMATION)
                    {
                        // Read file list
                        read_file_list();

                        // Read grid information
                        read_grid_header();

                        this->extent.resize(3);
                        this->extent[0].first = this->extent[1].first = this->extent[2].first = 0;
                        this->extent[0].second = static_cast<std::size_t>(this->header.x_resolution - 1);
                        this->extent[1].second = static_cast<std::size_t>(this->header.y_resolution - 1);
                        this->extent[2].second = static_cast<std::size_t>(this->header.z_resolution - 1);

                        // Read time steps and set number of components
                        this->timesteps.resize(this->data_filenames.size());

                        for (std::size_t i = 0; i < this->data_filenames.size(); ++i)
                        {
                            const auto data_header = read_data_header(i);
                            this->timesteps[i] = data_header.time_step;
                            this->num_components = static_cast<std::size_t>(data_header.components);
                            this->dataset_name = data_header.name + data_header.unit;
                        }

                        this->state = GLOBAL_INFORMATION;
                    }

                    if (this->state == EXTENT_INFORMATION)
                    {
                        // Read grid
                        read_grid();

                        this->node_coordinates.clear();
                        this->cell_sizes.clear();

                        this->node_coordinates.resize(3);
                        this->cell_sizes.resize(3);

                        for (std::size_t d = 0; d < 3; ++d)
                        {
                            // Calculate node coordinates
                            this->node_coordinates[d].reserve(this->cell_coordinates[d].size() + 1);
                            this->node_coordinates[d].push_back(static_cast<point_t>(0.0L));

                            for (std::size_t i = 1; i < this->cell_coordinates[d].size() + 1; ++i)
                            {
                                this->node_coordinates[d].push_back(static_cast<point_t>(2.0L) * this->cell_coordinates[d][i - 1] - this->node_coordinates[d][i - 1]);
                            }

                            // Calculate cell sizes
                            this->cell_sizes[d].reserve(this->cell_coordinates[d].size());

                            for (std::size_t i = 0; i < this->cell_coordinates[d].size(); ++i)
                            {
                                this->cell_sizes[d].push_back(this->node_coordinates[d][i + 1] - this->node_coordinates[d][i]);
                            }
                        }

                        // Extract subextent
                        for (std::size_t d = 0; d < 3; ++d)
                        {
                            if (this->subextent[d].first == this->extent[d].first)
                            {
                                this->cell_coordinates[d].resize(this->subextent[d].second - this->subextent[d].first + 1);
                                this->node_coordinates[d].resize(this->subextent[d].second - this->subextent[d].first + 2);
                                this->cell_sizes[d].resize(this->subextent[d].second - this->subextent[d].first + 1);
                            }
                            else
                            {
                                std::vector<point_t> sub_cell_coordinates, sub_node_coordinates, sub_cell_sizes;
                                sub_cell_coordinates.reserve(this->subextent[d].second - this->subextent[d].first + 1);
                                sub_node_coordinates.reserve(this->subextent[d].second - this->subextent[d].first + 2);
                                sub_cell_sizes.reserve(this->subextent[d].second - this->subextent[d].first + 1);

                                sub_cell_coordinates.insert(sub_cell_coordinates.end(), this->cell_coordinates[d].begin() + this->subextent[d].first,
                                    this->cell_coordinates[d].begin() + this->subextent[d].second + 1);
                                sub_node_coordinates.insert(sub_node_coordinates.end(), this->node_coordinates[d].begin() + this->subextent[d].first,
                                    this->node_coordinates[d].begin() + this->subextent[d].second + 2);
                                sub_cell_sizes.insert(sub_cell_sizes.end(), this->cell_sizes[d].begin() + this->subextent[d].first,
                                    this->cell_sizes[d].begin() + this->subextent[d].second + 1);

                                this->cell_coordinates[d].swap(sub_cell_coordinates);
                                this->node_coordinates[d].swap(sub_node_coordinates);
                                this->cell_sizes[d].swap(sub_cell_sizes);
                            }
                        }

                        this->state = ALL_INFORMATION;
                    }
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_reader<value_t, point_t>::read_file_list()
            {
                // Extract path from list file name
                const auto pos = list_filename.find_last_of("/\\");
                const auto path = list_filename.substr(0, pos + 1);

                // Read file name list
                std::ifstream ifs(list_filename, std::ifstream::in);
                std::string line;

                if (ifs.good())
                {
                    while (!ifs.eof())
                    {
                        std::getline(ifs, line);

                        if (!(line = boost::trim_copy(line)).empty())
                        {
                            this->data_filenames.push_back(path + line);
                        }
                    }

                    ifs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading file list!"));
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_reader<value_t, point_t>::read_grid_header()
            {
                // Create and read output header
                const std::size_t max_string = 80;

                std::ifstream ifs(this->grid_filename, std::ofstream::in | std::ofstream::binary);

                if (ifs.good())
                {
                    std::array<char, max_string> unit;

                    ifs.read(unit.data(), max_string);
                    ifs.read(reinterpret_cast<char*>(&this->header.x_resolution), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&this->header.y_resolution), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&this->header.z_resolution), sizeof(int));

                    this->header.unit.assign(unit.begin(), unit.end());
                    boost::trim(this->header.unit);

                    ifs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading grid file header!"));
                }
            }

            template <typename value_t, typename point_t>
            inline void fs3d_data_reader<value_t, point_t>::read_grid()
            {
                // Read data
                std::ifstream ifs(this->grid_filename, std::ifstream::in | std::ifstream::binary);

                // Create output
                this->cell_coordinates.clear();
                this->cell_coordinates.resize(3);

                this->cell_coordinates[0].resize(this->header.x_resolution);
                this->cell_coordinates[1].resize(this->header.y_resolution);
                this->cell_coordinates[2].resize(this->header.z_resolution);

                if (ifs.good())
                {
                    // Skip header
                    ifs.seekg(grid_header::length);

                    // Read for each dimension
                    if (std::is_same<typename grid_header::value_type, point_t>::value)
                    {
                        // Read block
                        ifs.read((char*)this->cell_coordinates[0].data(), this->cell_coordinates[0].size() * sizeof(typename grid_header::value_type));
                        ifs.read((char*)this->cell_coordinates[1].data(), this->cell_coordinates[1].size() * sizeof(typename grid_header::value_type));
                        ifs.read((char*)this->cell_coordinates[2].data(), this->cell_coordinates[2].size() * sizeof(typename grid_header::value_type));
                    }
                    else
                    {
                        for (std::size_t d = 0; d < 3; ++d)
                        {
                            // Read block to buffer
                            std::vector<typename grid_header::value_type> buffer(this->cell_coordinates[d].size());

                            ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(typename grid_header::value_type));

                            // Copy data from buffer
                            for (std::size_t i = 0; i < buffer.size(); ++i)
                            {
                                this->cell_coordinates[d][i] = static_cast<point_t>(buffer[i]);
                            }
                        }
                    }
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading grid file!"));
                }
            }

            template <typename value_t, typename point_t>
            inline data_header fs3d_data_reader<value_t, point_t>::read_data_header(const std::size_t timestep) const
            {
                // Create and read output header
                data_header header;
                const std::size_t max_string = 80;

                std::ifstream ifs(this->data_filenames[timestep], std::ofstream::in | std::ofstream::binary);

                if (ifs.good())
                {
                    std::array<char, max_string> name;
                    std::array<char, max_string> unit;

                    ifs.read(name.data(), max_string);
                    ifs.read(unit.data(), max_string);
                    ifs.read(reinterpret_cast<char*>(&header.time_id), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&header.time_step), sizeof(double));
                    ifs.read(reinterpret_cast<char*>(&header.x_resolution), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&header.y_resolution), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&header.z_resolution), sizeof(int));
                    ifs.read(reinterpret_cast<char*>(&header.identifier), sizeof(int));

                    header.name.assign(name.begin(), name.end());
                    boost::trim(header.name);

                    header.unit.assign(unit.begin(), unit.end());
                    boost::trim(header.unit);

                    ifs.seekg(0, std::ios::end);
                    const auto file_size = static_cast<std::size_t>(ifs.tellg());

                    const auto size_per_component = static_cast<std::size_t>(header.x_resolution) * static_cast<std::size_t>(header.y_resolution)
                        * static_cast<std::size_t>(header.z_resolution) * static_cast<std::size_t>(sizeof(double));

                    if (size_per_component == 0)
                    {
                        throw std::runtime_error(__tpf_error_message("Invalid resolution in header of file '", this->data_filenames[timestep], "'!"));
                    }

                    header.components = static_cast<int>((file_size - data_header::length) / size_per_component);

                    if (header.components * size_per_component != file_size - data_header::length)
                    {
                        throw std::runtime_error(__tpf_error_message("Invalid size of file '", this->data_filenames[timestep], "'!"));
                    }

                    ifs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading header of file '", this->data_filenames[timestep], "'!"));
                }

                return header;
            }

            template <typename value_t, typename point_t>
            template <std::size_t components>
            inline std::vector<value_t> fs3d_data_reader<value_t, point_t>::read_data(const std::size_t timestep) const
            {
                __tpf_timer_start(timer, "read_data", std::chrono::seconds);

                // Open file stream
                std::ifstream ifs(this->data_filenames[timestep], std::ifstream::in | std::ifstream::binary);

                // Calculate buffer and data size
                const std::size_t buffer_size = static_cast<std::size_t>(this->header.x_resolution) * 
                    static_cast<std::size_t>(this->header.y_resolution) * static_cast<std::size_t>(this->header.z_resolution);

                std::vector<value_t> data(buffer_size * components);

                if (ifs.good())
                {
                    // Skip header
                    ifs.seekg(data_header::length);

                    // Read for each component
                    if (std::is_same<typename data_header::value_type, value_t>::value && components == 1)
                    {
                        // Read block
                        ifs.read(reinterpret_cast<char*>(data.data()), data.size() * sizeof(typename data_header::value_type));
                    }
                    else
                    {
                        std::vector<typename data_header::value_type> buffer(buffer_size);

                        for (std::size_t c = 0; c < components; ++c)
                        {
                            // Read block
                            ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(typename data_header::value_type));

                            // Copy data from buffer
                            for (std::size_t i = 0; i < buffer_size; ++i)
                            {
                                const std::size_t data_index = (i * components) + c;

                                data[data_index] = static_cast<value_t>(buffer[i]);
                            }
                        }
                    }

                    ifs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading data file '", this->data_filenames[timestep], "'!"));
                }

                return data;
            }

            template <typename value_t, typename point_t>
            template <std::size_t components>
            inline std::vector<value_t> fs3d_data_reader<value_t, point_t>::read_data_block(const std::size_t timestep) const
            {
#ifdef __tpf_use_mpiio
                __tpf_timer("read_data_block parallel", std::chrono::seconds);

                const auto& mpi = mpi::get_instance();

                // Create sub array type representing the data block
                MPI_Datatype memtype, filetype;
                std::array<int, 4> memsizes, local_memsizes, local_memstarts;
                std::array<int, 3> filesizes, local_filesizes, local_filestarts;

                std::size_t memsize = 1;
                std::size_t filesize = 1;

                for (std::size_t i = 0; i < 3; ++i)
                {
                    memsizes[i] = utility::safe_downcast<int>(this->subextent[2 - i].second - this->subextent[2 - i].first + 1);
                    local_memsizes[i] = memsizes[i];
                    local_memstarts[i] = 0;
                    memsize *= memsizes[i];

                    filesizes[i] = utility::safe_downcast<int>(this->extent[2 - i].second - this->extent[2 - i].first + 1);
                    local_filesizes[i] = memsizes[i];
                    local_filestarts[i] = utility::safe_downcast<int>(this->subextent[2 - i].first - this->extent[2 - i].first);
                    filesize *= filesizes[i];
                }

                memsizes[3] = components;
                local_memsizes[3] = 1;
                local_memstarts[3] = 0;
                memsize *= components;

                MPI_Type_create_subarray(4, memsizes.data(), local_memsizes.data(), local_memstarts.data(), MPI_ORDER_C, MPI_DOUBLE, &memtype);
                MPI_Type_commit(&memtype);

                MPI_Type_create_subarray(3, filesizes.data(), local_filesizes.data(), local_filestarts.data(), MPI_ORDER_C, MPI_DOUBLE, &filetype);
                MPI_Type_commit(&filetype);

                // Open file
                std::vector<char> filename(this->data_filenames[timestep].begin(), this->data_filenames[timestep].end());
                filename.push_back('\0');
                MPI_File file;

                auto result = MPI_File_open(mpi.get_comm(), filename.data(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

                if (result != MPI_SUCCESS)
                {
                    throw std::runtime_error(__tpf_error_message("Unable to open data file '", this->data_filenames[timestep], "'!"));
                }

                // Calculate file size and offsets
                const MPI_Offset initial_offset = utility::safe_downcast<MPI_Offset>(data_header::length);
                const MPI_Offset increment = utility::safe_downcast<MPI_Offset>(filesize * sizeof(typename data_header::value_type));

                MPI_Offset file_size;

                MPI_File_get_size(file, &file_size);

                if (increment * components + initial_offset != file_size)
                {
                    throw std::runtime_error(__tpf_error_message("File size does not match for file '", this->data_filenames[timestep], "'!"));
                }

                // Read data block
                std::vector<typename data_header::value_type> data(memsize);

                for (std::size_t c = 0; c < components; ++c)
                {
                    result = MPI_File_set_view(file, initial_offset + c * increment, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

                    if (result != MPI_SUCCESS)
                    {
                        throw std::runtime_error(__tpf_error_message("Unable to set data file view for file '", this->data_filenames[timestep], "'!"));
                    }

                    result = MPI_File_read_all(file, &data[c], 1, memtype, MPI_STATUS_IGNORE);

                    if (result != MPI_SUCCESS)
                    {
                        throw std::runtime_error(__tpf_error_message("Unable to read data file for file '", this->data_filenames[timestep], "'!"));
                    }
                }

                MPI_Type_free(&memtype);
                MPI_Type_free(&filetype);

                MPI_File_close(&file);

                return utility::copy_data<value_t>(data);
#else
                __tpf_timer("read_data_block sequential", std::chrono::seconds);

                // Open file stream
                std::ifstream ifs(this->data_filenames[timestep], std::ifstream::in | std::ifstream::binary);

                // Create data array
                const std::array<std::size_t, 3> res = {
                    this->extent[0].second - this->extent[0].first + 1,
                    this->extent[1].second - this->extent[1].first + 1,
                    this->extent[2].second - this->extent[2].first + 1 };

                const std::array<std::size_t, 3> subres = {
                    this->subextent[0].second - this->subextent[0].first + 1,
                    this->subextent[1].second - this->subextent[1].first + 1,
                    this->subextent[2].second - this->subextent[2].first + 1 };

                std::vector<value_t> data(subres[0] * subres[1] * subres[2] * components);

                if (ifs.good())
                {
                    std::vector<typename data_header::value_type> buffer(subres[0]);

                    for (std::size_t c = 0; c < components; ++c)
                    {
                        for (std::size_t z = 0; z < subres[2]; ++z)
                        {
                            for (std::size_t y = 0; y < subres[1]; ++y)
                            {
                                // Skip to next block
                                const std::size_t file_index = this->subextent[0].first + res[0] *
                                    ((y + this->subextent[1].first) + res[1] * ((z + this->subextent[2].first) + res[2] * c));

                                ifs.seekg(data_header::length + file_index * sizeof(typename data_header::value_type));

                                if (std::is_same<typename data_header::value_type, value_t>::value && components == 1)
                                {
                                    // Read block
                                    const std::size_t data_index = c + components * (subres[0] * (y + subres[1] * z));

                                    ifs.read(reinterpret_cast<char*>(&data[data_index]), buffer.size() * sizeof(typename data_header::value_type));
                                }
                                else
                                {
                                    // Read block
                                    ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size() * sizeof(typename data_header::value_type));

                                    // Copy data from buffer
                                    for (std::size_t i = 0; i < subres[0]; ++i)
                                    {
                                        const std::size_t data_index = c + components * (i + subres[0] * (y + subres[1] * z));

                                        data[data_index] = static_cast<value_t>(buffer[i]);
                                    }
                                }
                            }
                        }
                    }

                    ifs.close();
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Error reading data file '", this->data_filenames[timestep], "'!"));
                }

                return data;
#endif
            }
        }
    }
}
