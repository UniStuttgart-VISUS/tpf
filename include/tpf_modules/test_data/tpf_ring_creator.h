#pragma once

#include "tpf_geometry_creator.h"

#include "tpf/data/tpf_grid.h"

#include "tpf/math/tpf_vector.h"

#include <cstddef>

namespace tpf
{
    namespace modules
    {
        namespace test_data_aux
        {
            /// <summary>
            /// Class for the creation of rings
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class ring_creator : public geometry_creator<floatp_t>
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="fractions">Fractions (VOF)</param>
                /// <param name="velocities">Velocities</param>
                /// <param name="num_cells">Number of grid cells</param>
                ring_creator(data::grid<floatp_t, floatp_t, 3, 1>* fractions, data::grid<floatp_t, floatp_t, 3, 3>* velocities, std::size_t num_cells);

                /// <summary>
                /// Add spinning velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of spinning</param>
                virtual void add_spinning(floatp_t magnitude = static_cast<floatp_t>(1.0L)) override;

                /// <summary>
                /// Add movement velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of movement</param>
                virtual void add_movement(floatp_t magnitude = static_cast<floatp_t>(1.0L)) override;

                /// <summary>
                /// Add expanding velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of expansion</param>
                virtual void add_expansion(floatp_t magnitude = static_cast<floatp_t>(1.0L)) override;

                /// <summary>
                /// Add rotating velocities
                /// </summary>
                /// <param name="magnitude">Magnitude of rotation</param>
                virtual void add_rotation(floatp_t magnitude = static_cast<floatp_t>(1.0L)) override;

            protected:
                /// <summary>
                /// Create a static orb
                /// </summary>
                virtual void create() override;

            private:
                /// Ring information
                const math::vec3_t<floatp_t> ring_center;
                const math::vec3_t<floatp_t> ring_plane_normal;

                const floatp_t ring_inner_radius;
                const floatp_t ring_outer_radius;
                const floatp_t ring_radius;
            };
        }
    }
}

#include "tpf_ring_creator.inl"
