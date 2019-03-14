#pragma once

#include "tpf_geometric_float.h"
#include "tpf_plane.h"
#include "tpf_point.h"

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Calculate distance between two points
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="first_point">First point</param>
        /// <param name="second_point">Second point</param>
        /// <returns>Distance</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static geometric_float<floatp_t, typename kernel_t::FT> calculate_squared_distance(const point<floatp_t, kernel_t>& first_point, const point<floatp_t, kernel_t>& second_point);

        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static floatp_t calculate_distance(const point<floatp_t, kernel_t>& first_point, const point<floatp_t, kernel_t>& second_point);

        /// <summary>
        /// Calculate the distance between a point and a plane
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="point">Point</param>
        /// <param name="plane">Plane</param>
        /// <returns>Distance</returns>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static geometric_float<floatp_t, typename kernel_t::FT> calculate_squared_distance(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane);

        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static floatp_t calculate_distance(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane);
    }
}

#include "tpf_distance.inl"