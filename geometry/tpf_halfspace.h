#pragma once

#include "tpf_plane.h"
#include "tpf_point.h"

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Check for point position relative to plane
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="plane">Plane</param>
        /// <returns>True: point is on the plane</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool is_on_plane(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane);

        /// <summary>
        /// Check for point position relative to plane
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="plane">Plane</param>
        /// <returns>True: point is in positive halfspace, otherwise on the plane or in negative halfspace</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool is_in_positive_halfspace(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane);

        /// <summary>
        /// Check for point position relative to plane
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="plane">Plane</param>
        /// <returns>True: point is in positive halfspace, otherwise on the plane or in negative halfspace</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool is_in_negative_halfspace(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane);
    }
}

#include "tpf_halfspace.inl"