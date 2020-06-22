#pragma once

#include "tpf_fs3d_data.h"

#include "tpf/module/tpf_module_reader_base.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include <string>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_reader_aux
        {
            /// <summary>
            /// Class to provide data from  files
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            template <typename value_t, typename point_t>
            class fs3d_data_reader
            {
                static_assert(sizeof(double) == 8, "double needs to be of 8 bytes (64 bits).");
                static_assert(sizeof(int) == 4, "int needs to be of 4 bytes (32 bits).");

            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="list_filename">List file name</param>
                /// <param name="grid_filename">Grid file name</param>
                fs3d_data_reader(const std::string& list_filename, const std::string& grid_filename);

                /// <summary>
                /// Return information on spatial and temporal extents
                /// </summary>
                /// <returns>Spatial and temporal extents</returns>
                data_information get_information();

                /// <summary>
                /// Set sub extent
                /// </summary>
                /// <param name="extent">Local subextent</param>
                void set_subextent(const data::extent_t& extent);

                /// <summary>
                /// Read one time step and return corresponding grid
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <returns>Loaded grid</returns>
                /// <template name="components">Number of components</template>
                template <std::size_t components>
                data::grid<value_t, point_t, 3, components> read_timestep(std::size_t timestep);

            private:
                /// <summary>
                /// Read information from grid and data files
                /// </summary>
                void read_information();

                /// <summary>
                /// Read list file containing all data file names
                /// </summary>
                void read_file_list();

                /// <summary>
                /// Read grid file header
                /// </summary>
                void read_grid_header();

                /// <summary>
                /// Read grid file
                /// </summary>
                void read_grid();

                /// <summary>
                /// Read data file header
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <returns>Data file header</returns>
                data_header read_data_header(std::size_t timestep) const;

                /// <summary>
                /// Read data file
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <returns>Data</returns>
                /// <template name="components">Number of components</template>
                template <std::size_t components>
                std::vector<value_t> read_data(std::size_t timestep) const;

                /// <summary>
                /// Read data block
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <returns>Data block</returns>
                /// <template name="components">Number of components</template>
                template <std::size_t components>
                std::vector<value_t> read_data_block(std::size_t timestep) const;

                /// List and grid file name
                const std::string list_filename;
                const std::string grid_filename;

                /// Data file names
                std::vector<std::string> data_filenames;

                /// Time steps available
                std::vector<double> timesteps;

                /// Grid
                grid_header header;

                /// Implicit grid information
                data::extent_t extent;
                data::extent_t subextent;

                typename data::grid_information<point_t>::array_type cell_coordinates;
                typename data::grid_information<point_t>::array_type node_coordinates;
                typename data::grid_information<point_t>::array_type cell_sizes;

                /// Data set name
                std::string dataset_name;

                /// Number of components
                std::size_t num_components;

                /// State of information
                enum state_t
                {
                    NO_INFORMATION, GLOBAL_INFORMATION, EXTENT_INFORMATION, ALL_INFORMATION
                } state;
            };
        }
    }
}

#include "tpf_fs3d_data_reader.inl"