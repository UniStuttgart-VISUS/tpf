#pragma once

#include <string>

namespace tpf
{
    namespace modules
    {
        namespace fs3d_reader_aux
        {
            /// <summary>
            /// Header for data files
            /// </summary>
            struct data_header
            {
                /// Header length in bytes
                static const std::size_t length = 188;

                /// Array name
                std::basic_string<char> name;

                /// Unit of stored values
                std::basic_string<char> unit;

                /// Time step information
                int time_id;
                double time_step;

                /// Grid resolution
                int x_resolution;
                int y_resolution;
                int z_resolution;

                /// Identifier (usually 1)
                int identifier;

                /// Number of components (implicit from file size)
                int components;

                /// Value type of data read from data files
                typedef double value_type;
            };

            /// <summary>
            /// Header for grid files
            /// </summary>
            struct grid_header
            {
                /// Header length in bytes
                static const std::size_t length = 92;

                /// Unit of stored values
                std::basic_string<char> unit;

                /// Grid resolution
                int x_resolution;
                int y_resolution;
                int z_resolution;

                /// Value type of data read from grid files
                typedef double value_type;
            };
        }
    }
}