#include "tpf_rectangle.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../log/tpf_log.h"

#include <CGAL/Point_2.h>
#include <CGAL/Iso_rectangle_2.h>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>::rectangle(const point<floatp_t, kernel_t>& point_min, const point<floatp_t, kernel_t>& point_max, const bool render_as_triangles) noexcept
        {
            const typename kernel_t::Point_2 bottom_left_point(point_min.get_vertex()[0], point_min.get_vertex()[1]);
            const typename kernel_t::Point_2 top_right_point(point_max.get_vertex()[0], point_max.get_vertex()[1]);

            this->_rectangle = typename kernel_t::Iso_rectangle_2(bottom_left_point, top_right_point);

            this->render_as_triangles = render_as_triangles;
        }

        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>::rectangle(const rectangle<floatp_t, kernel_t>& copy) noexcept
        {
            this->_rectangle = copy._rectangle;

            this->render_as_triangles = copy.render_as_triangles;
        }

        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>::rectangle(const typename kernel_t::Iso_rectangle_2& copy) noexcept
        {
            this->_rectangle = copy;

            this->render_as_triangles = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>& rectangle<floatp_t, kernel_t>::operator=(const rectangle<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_rectangle = copy._rectangle;

            this->render_as_triangles = copy.render_as_triangles;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>& rectangle<floatp_t, kernel_t>::operator=(const typename kernel_t::Iso_rectangle_2& copy) noexcept
        {
            this->_rectangle = copy;

            this->render_as_triangles = false;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool rectangle<floatp_t, kernel_t>::operator==(const rectangle<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_rectangle == other._rectangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool rectangle<floatp_t, kernel_t>::operator==(const typename kernel_t::Iso_rectangle_2& other) const noexcept
        {
            return this->_rectangle == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> rectangle<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<rectangle<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& rectangle<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                throw exception::not_implemented_exception(__tpf_error_message("Transformation not available for 2D objects."));
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> rectangle<floatp_t, kernel_t>::calculate_area() const
        {
            return this->_rectangle.area();
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> rectangle<floatp_t, kernel_t>::get_min_point() const
        {
            return this->_rectangle.min();
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> rectangle<floatp_t, kernel_t>::get_max_point() const
        {
            return this->_rectangle.max();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t rectangle<floatp_t, kernel_t>::get_num_points() const
        {
            return 4;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> rectangle<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;

            for (std::size_t i = 0; i < 4; ++i)
            {
                const auto point_2d = this->_rectangle.vertex(static_cast<int>(i));

                points.emplace_back(static_cast<floatp_t>(CGAL::to_double(point_2d[0])),
                    static_cast<floatp_t>(CGAL::to_double(point_2d[1])), static_cast<floatp_t>(0.0));
            }

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t rectangle<floatp_t, kernel_t>::get_num_cells() const
        {
            if (this->render_as_triangles)
            {
                return 2;
            }
            else
            {
                return 1;
            }
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> rectangle<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells;

            if (this->render_as_triangles)
            {
                cells.resize(2);
                cells[0].resize(3);
                cells[1].resize(3);

                cells[0][0] = 0; cells[0][1] = 1; cells[0][2] = 2;
                cells[1][0] = 2; cells[1][1] = 3; cells[1][2] = 0;
            }
            else
            {
                cells.resize(1);
                cells[0].resize(4);
                
                std::iota(cells[0].begin(), cells[0].end(), static_cast<std::size_t>(0));
            }

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> rectangle<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> rectangle<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Iso_rectangle_2& rectangle<floatp_t, kernel_t>::get_internal() const
        {
            return this->_rectangle;
        }

        template <typename floatp_t, typename kernel_t>
        inline rectangle<floatp_t, kernel_t>::operator const typename kernel_t::Iso_rectangle_2&() const
        {
            return this->_rectangle;
        }
    }
}
