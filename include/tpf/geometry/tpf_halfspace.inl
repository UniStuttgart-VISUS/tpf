#include "tpf_halfspace.h"

#include "tpf_plane.h"
#include "tpf_point.h"

#include <CGAL/enum.h>
#include <CGAL/squared_distance_3.h>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline bool is_on_plane(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane)
        {
            return plane.get_internal().has_on(point.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline bool is_in_positive_halfspace(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane)
        {
            return plane.get_internal().oriented_side(point.get_internal()) == CGAL::Oriented_side::ON_POSITIVE_SIDE;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool is_in_negative_halfspace(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane)
        {
            return plane.get_internal().oriented_side(point.get_internal()) == CGAL::Oriented_side::ON_NEGATIVE_SIDE;
        }
    }
}