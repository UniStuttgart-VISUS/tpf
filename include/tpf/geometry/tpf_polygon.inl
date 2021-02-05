#include "tpf_polygon.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_triangle.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../log/tpf_log.h"

#include "../math/tpf_transformer.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::polygon(const std::vector<point<floatp_t, kernel_t>>& points, const bool make_convex, const bool render_as_triangles) noexcept(false)
            : polygon(points, calculate_normal(points), make_convex, render_as_triangles)
        { }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::polygon(std::vector<point<floatp_t, kernel_t>> points, const Eigen::Matrix<floatp_t, 3, 1>& normal,
            const bool make_convex, const bool render_as_triangles) noexcept(false) : render_as_triangles(render_as_triangles)
        {
            if (points.size() < 3)
            {
                throw std::runtime_error(__tpf_error_message("A polygon cannot be constructed from less than 3 points."));
            }

            // Transform points to 2D
            this->transformation = tpf::math::transformer<floatp_t, 3>(points[0].get_vertex(), normal);

            std::vector<typename kernel_t::Point_2> points_2D;
            
            for (const auto& point : points)
            {
                const auto point_2D = this->transformation.transform_inverse(point.get_vertex());

                points_2D.push_back(typename kernel_t::Point_2(point_2D[0], point_2D[1]));
            }

            // Calculate convex hull
            if (make_convex)
            {
                std::vector<typename kernel_t::Point_2> convex_hull;
                CGAL::convex_hull_2(points_2D.begin(), points_2D.end(), std::back_inserter(convex_hull));

                std::swap(points_2D, convex_hull);
            }

            // Save polygon
            this->_polygon = CGAL::Polygon_2<kernel_t>(points_2D.begin(), points_2D.end());
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::polygon(const triangle<floatp_t, kernel_t>& triangle) noexcept
            : polygon({ triangle.get_internal()[0], triangle.get_internal()[1], triangle.get_internal()[2] }, false, true)
        { }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::polygon(const CGAL::Polygon_2<kernel_t>& polygon) noexcept
        {
            this->_polygon = polygon;

            this->render_as_triangles = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::polygon(const polygon<floatp_t, kernel_t>& copy) noexcept
        {
            this->_polygon = copy._polygon;

            this->transformation = copy.transformation;

            this->render_as_triangles = copy.render_as_triangles;
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>& polygon<floatp_t, kernel_t>::operator=(const triangle<floatp_t, kernel_t>& triangle) noexcept
        {
            return operator=(polygon<floatp_t, kernel_t>({ triangle.get_internal()[0], triangle.get_internal()[1], triangle.get_internal()[2] }, false, true));
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>& polygon<floatp_t, kernel_t>::operator=(const CGAL::Polygon_2<kernel_t>& polygon) noexcept
        {
            this->_polygon = polygon;

            this->transformation = math::transformer<floatp_t, 3>();

            this->render_as_triangles = false;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>& polygon<floatp_t, kernel_t>::operator=(const polygon<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_polygon = copy._polygon;

            this->transformation = copy.transformation;

            this->render_as_triangles = copy.render_as_triangles;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool polygon<floatp_t, kernel_t>::operator==(const polygon<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_polygon == other._polygon && this->transformation == other.transformation;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool polygon<floatp_t, kernel_t>::operator==(const CGAL::Polygon_2<kernel_t>& other) const noexcept
        {
            return this->_polygon == other && this->transformation.is_unit();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> polygon<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<polygon<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& polygon<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                throw exception::not_implemented_exception(__tpf_error_message("Transformation not available for 2D objects."));
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> polygon<floatp_t, kernel_t>::calculate_area() const
        {
            return this->_polygon.area();
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> polygon<floatp_t, kernel_t>::calculate_centroid() const
        {
            // Check convexity
            if (!is_convex())
            {
                throw exception::not_implemented_exception(__tpf_error_message("Calculation of centroid only implemented for convex polygons."));
            }

            // Get transformed points
            const auto points = get_points();

            // Calculate centroid
            typename kernel_t::FT temp, A, c_x, c_y;
            A = c_x = c_y = typename kernel_t::FT(0.0);

            for (auto it = this->_polygon.vertices_begin(); it != this->_polygon.vertices_end(); ++it)
            {
                const auto it_next = (std::next(it) == this->_polygon.vertices_end()) ? this->_polygon.vertices_begin() : std::next(it);

                temp = (*it)[0] * (*it_next)[1] - (*it_next)[0] * (*it)[1];

                c_x += ((*it)[0] + (*it_next)[0]) * temp;
                c_y += ((*it)[1] + (*it_next)[1]) * temp;

                A += temp;
            }

            c_x /= typename kernel_t::FT(3.0) * A;
            c_y /= typename kernel_t::FT(3.0) * A;

            return this->transformation.transform(Eigen::Matrix<floatp_t, 3, 1>(static_cast<floatp_t>(CGAL::to_double(c_x)),
                static_cast<floatp_t>(CGAL::to_double(c_y)), static_cast<floatp_t>(0.0)));
        }

        template <typename floatp_t, typename kernel_t>
        bool polygon<floatp_t, kernel_t>::is_convex() const
        {
            return this->_polygon.is_convex();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t polygon<floatp_t, kernel_t>::get_num_points() const
        {
            return this->_polygon.size();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> polygon<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.reserve(get_num_points());

            for (auto it = this->_polygon.vertices_begin(); it != this->_polygon.vertices_end(); ++it)
            {
                points.push_back(this->transformation.transform(Eigen::Matrix<floatp_t, 3, 1>(
                    static_cast<floatp_t>(CGAL::to_double((*it)[0])), static_cast<floatp_t>(CGAL::to_double((*it)[1])), static_cast<floatp_t>(0.0))));
            }

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t polygon<floatp_t, kernel_t>::get_num_cells() const
        {
            if (this->render_as_triangles && this->_polygon.is_convex())
            {
                return this->_polygon.size() - 2;
            }

            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> polygon<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells;

            if (this->render_as_triangles && this->_polygon.is_convex())
            {
                cells.resize(get_num_cells());

                for (std::size_t c = 0; c < cells.size(); ++c)
                {
                    cells[c].push_back(0);
                    cells[c].push_back(c + 1);
                    cells[c].push_back(c + 2);
                }
            }
            else
            {
                cells.resize(static_cast<std::size_t>(1));

                for (std::size_t i = 0; i < get_num_points(); ++i)
                {
                    cells.front().push_back(i);
                }
            }

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> polygon<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> polygon<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const CGAL::Polygon_2<kernel_t>& polygon<floatp_t, kernel_t>::get_internal() const
        {
            return this->_polygon;
        }

        template <typename floatp_t, typename kernel_t>
        inline polygon<floatp_t, kernel_t>::operator const CGAL::Polygon_2<kernel_t>&() const
        {
            return this->_polygon;
        }

        template <typename floatp_t, typename kernel_t>
        inline const math::transformer<floatp_t, 3>& polygon<floatp_t, kernel_t>::get_transformation() const
        {
            return this->transformation;
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> polygon<floatp_t, kernel_t>::calculate_normal(const std::vector<point<floatp_t, kernel_t>> &points) const noexcept(false)
        {
            if (points.size() < 3)
            {
                throw std::runtime_error(__tpf_error_message("A polygon cannot be constructed from less than 3 points."));
            }

            // Calculate normal
            const auto AB = points[1].get_internal() - points[0].get_internal();
            const auto AC = points[2].get_internal() - points[0].get_internal();

            const auto direction = CGAL::cross_product(AB, AC);

            return point<floatp_t, kernel_t>(CGAL::ORIGIN + (direction / std::sqrt(CGAL::to_double(direction.squared_length()))));
        }
    }
}
