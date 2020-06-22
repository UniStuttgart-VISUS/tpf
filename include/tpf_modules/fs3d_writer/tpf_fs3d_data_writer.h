#pragma once

#include "../fs3d_reader/tpf_fs3d_data.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include <string>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_writer_aux
        {
            /// <summary>
            /// Class to provide data from  files
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            template <typename value_t, typename point_t>
            class fs3d_data_writer
            {
                static_assert(sizeof(double) == 8, "double needs to be of 8 bytes (64 bits).");
                static_assert(sizeof(int) == 4, "int needs to be of 4 bytes (32 bits).");

            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="data_filename">Data file name</param>
                /// <param name="grid_filename">Grid file name</param>
                fs3d_data_writer(const std::string& data_filename, const std::string& grid_filename);

                /// <summary>
                /// Set global extent
                /// </summary>
                /// <param name="global_extent">Global extent</param>
                void set_global_extent(const data::extent_t& global_extent);

                /// <summary>
                /// Write one time step
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <param name="time">Time</param>
                /// <param name="grid">Grid to write</param>
                /// <template name="Components">Number of components</template>
                template <int Components>
                void write_timestep(std::size_t timestep, double time, const std::string& name, const std::string& unit,
                    const std::string& grid_unit, const data::grid<value_t, point_t, 3, Components>& grid);

            private:
                /// <summary>
                /// Write grid file header
                /// </summary>
                void write_grid_header(const fs3d_reader_aux::grid_header& header) const;

                /// <summary>
                /// Write grid file
                /// </summary>
                void write_grid(const data::grid_information<point_t>& grid) const;

                /// <summary>
                /// Write data file header
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <param name="header">Data file header</param>
                void write_data_header(std::size_t timestep, const fs3d_reader_aux::data_header& header) const;

                /// <summary>
                /// Write data file
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <param name="data">Data</param>
                /// <template name="Components">Number of components</template>
                template <int Components>
                void write_data(std::size_t timestep, const std::vector<value_t>& data) const;

                /// <summary>
                /// Write data block
                /// </summary>
                /// <param name="timestep">Time step</param>
                /// <param name="data">Data block</data>
                /// <template name="Components">Number of components</template>
                template <int Components>
                void write_data_block(std::size_t timestep, const std::vector<value_t>& data) const;

                /// Grid file name
                const std::string grid_filename;
                const std::string data_filename;

                /// Implicit grid information
                data::extent_t global_extent;
                data::extent_t subextent;
            };
        }
    }
}

#include "tpf_fs3d_data_writer.inl"
