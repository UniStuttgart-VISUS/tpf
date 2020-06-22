#pragma once

#include "tpf/data/tpf_grid.h"

#include <cstddef>

namespace tpf
{
    namespace modules
    {
        namespace test_data_aux
        {
            /// <summary>
            /// Base class for the creation of geometrical objects
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class geometry_creator
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="fractions">Fractions (VOF)</param>
                /// <param name="velocities">Velocities</param>
                /// <param name="num_cells">Number of grid cells</param>
                geometry_creator(data::grid<floatp_t, floatp_t, 3, 1>* fractions, data::grid<floatp_t, floatp_t, 3, 3>* velocities, std::size_t num_cells);

                /// <summary>
                /// Add spinning velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of spinning</param>
                virtual void add_spinning(floatp_t magnitude) = 0;

                /// <summary>
                /// Add movement velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of movement</param>
                virtual void add_movement(floatp_t magnitude) = 0;

                /// <summary>
                /// Add expanding velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of expansion</param>
                virtual void add_expansion(floatp_t magnitude) = 0;

                /// <summary>
                /// Add rotating velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of rotation</param>
                virtual void add_rotation(floatp_t magnitude) = 0;

            protected:
                /// <summary>
                /// Create a static object
                /// </summary>
                virtual void create() = 0;

                /// Grid information
                const std::size_t num_cells;

                /// Fractions and velocities
                data::grid<floatp_t, floatp_t, 3, 1>* fractions;
                data::grid<floatp_t, floatp_t, 3, 3>* velocities;

                /// Store positions temporarily
                data::grid<floatp_t, floatp_t, 3, 3> positions;
            };
        }
    }
}

#include "tpf_geometry_creator.inl"
