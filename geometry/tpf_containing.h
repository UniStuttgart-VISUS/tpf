#pragma once

#include "tpf_cuboid.h"
#include "tpf_point.h"
#include "tpf_rectangle.h"

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Check for point position relative to rectangle
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="rectangle">Rectangle</param>
        /// <returns>True: point is inside the rectangle</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool is_in_rectangle(const point<floatp_t, kernel_t>& point, const rectangle<floatp_t, kernel_t>& rectangle);

        /// <summary>
        /// Check for point position relative to cuboid
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="cuboid">Cuboid</param>
        /// <returns>True: point is inside the cuboid</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool is_in_cuboid(const point<floatp_t, kernel_t>& point, const cuboid<floatp_t, kernel_t>& cuboid);
    }
}

#include "tpf_containing.inl"