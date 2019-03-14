#include "tpf_triangle.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../exception/tpf_not_implemented_exception.h"

#include <CGAL/Triangle_3.h>

#include <cmath>
#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>::triangle(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2, const point<floatp_t, kernel_t>& point_3) noexcept
        {
            this->_triangle = typename kernel_t::Triangle_3(point_1.get_internal(), point_2.get_internal(), point_3.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>::triangle(const typename kernel_t::Triangle_3& triangle) noexcept
        {
            this->_triangle = triangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>::triangle(const triangle<floatp_t, kernel_t>& copy) noexcept
        {
            this->_triangle = copy._triangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>& triangle<floatp_t, kernel_t>::operator=(const typename kernel_t::Triangle_3& triangle) noexcept
        {
            geometric_object<floatp_t>::operator=(triangle);

            this->_triangle = triangle;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>& triangle<floatp_t, kernel_t>::operator=(const triangle<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_triangle = copy._triangle;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool triangle<floatp_t, kernel_t>::operator==(const triangle<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_triangle == other._triangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool triangle<floatp_t, kernel_t>::operator==(const typename kernel_t::Triangle_3& other) const noexcept
        {
            return this->_triangle == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> triangle<floatp_t, kernel_t>::clone() const
        {
            return std::make_shared<triangle<floatp_t, kernel_t>>(*this);
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> triangle<floatp_t, kernel_t>::calculate_squared_area() const
        {
            return this->_triangle.squared_area();
        }

        template <typename floatp_t, typename kernel_t>
        inline floatp_t triangle<floatp_t, kernel_t>::calculate_area() const
        {
            return static_cast<floatp_t>(std::sqrt(CGAL::to_double(this->_triangle.squared_area())));
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t triangle<floatp_t, kernel_t>::get_num_points() const
        {
            return 3;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> triangle<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;

            points.push_back(point<floatp_t, kernel_t>(this->_triangle.vertex(0)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_triangle.vertex(1)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_triangle.vertex(2)).get_vertex());
            
            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t triangle<floatp_t, kernel_t>::get_num_cells() const
        {
            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> triangle<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(1));
            cells.front().push_back(static_cast<std::size_t>(0));
            cells.front().push_back(static_cast<std::size_t>(1));
            cells.front().push_back(static_cast<std::size_t>(2));

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> triangle<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> triangle<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Triangle_3& triangle<floatp_t, kernel_t>::get_internal() const
        {
            return this->_triangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline triangle<floatp_t, kernel_t>::operator const typename kernel_t::Triangle_3&() const
        {
            return this->_triangle;
        }
    }
}
