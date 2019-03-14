#pragma once

#include "tpf_cuboid.h"
#include "tpf_line.h"
#include "tpf_plane.h"
#include "tpf_point.h"
#include "tpf_polyhedron.h"
#include "tpf_tetrahedron.h"

#include "../utility/tpf_optional.h"

#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Does the line intersect with a plane?
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="line">Line</param>
        /// <param name="plane">Plane</param>
        /// <return>Is there an intersection?</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool does_intersect_with(const line<floatp_t, kernel_t>& line, const plane<floatp_t, kernel_t>& plane);

        /// <summary>
        /// Does the plane intersect with a cuboid?
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="plane">Plane</param>
        /// <param name="cuboid">Cuboid</param>
        /// <return>Is there an intersection?</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static bool does_intersect_with(const plane<floatp_t, kernel_t>& plane, const cuboid<floatp_t, kernel_t>& cuboid);

        // ------------------------------------------------------------------------------------------------------------------------------

        /// <summary>
        /// Intersect the line with a plane
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="line">Line</param>
        /// <param name="plane">Plane</param>
        /// <return>Intersection</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static utility::optional<point<floatp_t, kernel_t>> intersect_with(const line<floatp_t, kernel_t>& line, const plane<floatp_t, kernel_t>& plane);

        /// <summary>
        /// Intersect the plane with a cuboid
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="plane">Plane</param>
        /// <param name="cuboid">Cuboid</param>
        /// <return>Intersections</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const cuboid<floatp_t, kernel_t>& cuboid);

        /// <summary>
        /// Intersect the plane with a tetrahedron
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="plane">Plane</param>
        /// <param name="tetrahedron">Tetrahedron</param>
        /// <return>Intersections</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const tetrahedron<floatp_t, kernel_t>& tetrahedron);

        /// <summary>
        /// Intersect the plane with a polyhedron
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">kernel_t type</template>
        /// <param name="plane">Plane</param>
        /// <param name="polyhedron">Polyhedron</param>
        /// <return>Intersections</return>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const polyhedron<floatp_t, kernel_t>& polyhedron);
    }
}

#include "tpf_intersection.inl"