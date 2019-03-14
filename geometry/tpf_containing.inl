#include "tpf_containing.h"

#include "tpf_cuboid.h"
#include "tpf_point.h"
#include "tpf_rectangle.h"

#include <CGAL/enum.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        static bool is_in_rectangle(const point<floatp_t, kernel_t>& point, const rectangle<floatp_t, kernel_t>& rectangle)
        {
            const auto min_point_3 = rectangle.get_min_point().get_internal();
            const auto max_point_3 = rectangle.get_max_point().get_internal();

            const typename kernel_t::Point_2 min_point_2(min_point_3.x(), min_point_3.y());
            const typename kernel_t::Point_2 max_point_2(max_point_3.x(), max_point_3.y());

            const typename kernel_t::Point_2 point_2(point.get_internal().x(), point.get_internal().y());
            const typename kernel_t::Iso_rectangle_2 cgal_rectangle(min_point_2, max_point_2);

            return cgal_rectangle.bounded_side(point_2) != CGAL::Bounded_side::ON_UNBOUNDED_SIDE;
        }

        template <typename floatp_t, typename kernel_t>
        static bool is_in_cuboid(const point<floatp_t, kernel_t>& point, const cuboid<floatp_t, kernel_t>& cuboid)
        {
            return cuboid.get_internal().bounded_side(point.get_internal()) != CGAL::Bounded_side::ON_UNBOUNDED_SIDE;
        }
    }
}