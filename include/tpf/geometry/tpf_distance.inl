#include "tpf_distance.h"

#include "tpf_plane.h"
#include "tpf_point.h"

#include <CGAL/squared_distance_3.h>

#include <cmath>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> calculate_squared_distance(const point<floatp_t, kernel_t>& first_point, const point<floatp_t, kernel_t>& second_point)
        {
            return CGAL::squared_distance(first_point.get_internal(), second_point.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline floatp_t calculate_distance(const point<floatp_t, kernel_t>& first_point, const point<floatp_t, kernel_t>& second_point)
        {
            return std::sqrt(calculate_squared_distance(first_point, second_point).get_float_value());
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> calculate_squared_distance(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane)
        {
            return CGAL::squared_distance(plane.get_internal(), point.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline floatp_t calculate_distance(const point<floatp_t, kernel_t>& point, const plane<floatp_t, kernel_t>& plane)
        {
            return std::sqrt(calculate_squared_distance(point, plane).get_float_value());
        }
    }
}