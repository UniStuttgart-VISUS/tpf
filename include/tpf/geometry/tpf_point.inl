#include "tpf_point.h"

#include "tpf_geometric_object.h"

#include "../log/tpf_log.h"

#include <CGAL/Point_3.h>

#include <iterator>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::point(floatp_t x, floatp_t y, floatp_t z) noexcept
        {
            this->_point = typename kernel_t::Point_3(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z));
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::point(const Eigen::Matrix<floatp_t, 3, 1>& point) noexcept
        {
            this->_point = typename kernel_t::Point_3(static_cast<double>(point[0]), static_cast<double>(point[1]), static_cast<double>(point[2]));
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::point(const typename kernel_t::Point_3& point) noexcept
        {
            this->_point = point;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::point(const typename kernel_t::Point_2& point, floatp_t z) noexcept
        {
            this->_point = typename kernel_t::Point_3(point[0], point[1], z);
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::point(const point<floatp_t, kernel_t>& copy) noexcept
        {
            this->_point = copy._point;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>& point<floatp_t, kernel_t>::operator=(const typename kernel_t::Point_3& point) noexcept
        {
            geometric_object<floatp_t>::operator=(point);

            this->_point = point;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>& point<floatp_t, kernel_t>::operator=(const typename kernel_t::Point_2& point) noexcept
        {
            geometric_object<floatp_t>::operator=(point);

            this->_point = typename kernel_t::Point_3(point[0], point[1], 0.0);

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>& point<floatp_t, kernel_t>::operator=(const point<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_point = copy._point;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>& point<floatp_t, kernel_t>::operator=(const Eigen::Matrix<floatp_t, 3, 1>& point) noexcept
        {
            this->_point = typename kernel_t::Point_3(static_cast<double>(point[0]), static_cast<double>(point[1]), static_cast<double>(point[2]));

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool point<floatp_t, kernel_t>::operator==(const point<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_point == other._point;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool point<floatp_t, kernel_t>::operator==(const typename kernel_t::Point_3& other) const noexcept
        {
            return this->_point == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> point<floatp_t, kernel_t>::clone() const
        {
            return std::make_shared<point<floatp_t, kernel_t>>(*this);
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> point<floatp_t, kernel_t>::get_vertex() const
        {
            return Eigen::Matrix<floatp_t, 3, 1>(
                static_cast<floatp_t>(CGAL::to_double(this->_point.x())),
                static_cast<floatp_t>(CGAL::to_double(this->_point.y())),
                static_cast<floatp_t>(CGAL::to_double(this->_point.z())));
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::operator Eigen::Matrix<floatp_t, 3, 1>() const
        {
            return get_vertex();
        }

        template <typename floatp_t, typename kernel_t>
        inline floatp_t point<floatp_t, kernel_t>::operator[](std::size_t i) const
        {
#ifdef __tpf_sanity_checks
            if (i > 2)
            {
                throw std::runtime_error(__tpf_error_message("Illegal vector element. Index must be 0, 1, or 2."));
            }
#endif

            return static_cast<floatp_t>(CGAL::to_double(this->_point[static_cast<int>(i)]));
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t point<floatp_t, kernel_t>::get_num_points() const
        {
            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> point<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.push_back(get_vertex());

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t point<floatp_t, kernel_t>::get_num_cells() const
        {
            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> point<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(1));
            cells.front().push_back(static_cast<std::size_t>(0));

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> point<floatp_t, kernel_t>::serialize() const
        {
            std::vector<char> buffer(1 + 3 * sizeof(floatp_t));

            buffer[0] = static_cast<char>(geometry_t::POINT);

            const auto vertex = get_vertex();

            *reinterpret_cast<floatp_t*>(&buffer[1 + 0 * sizeof(floatp_t)]) = vertex[0];
            *reinterpret_cast<floatp_t*>(&buffer[1 + 1 * sizeof(floatp_t)]) = vertex[1];
            *reinterpret_cast<floatp_t*>(&buffer[1 + 2 * sizeof(floatp_t)]) = vertex[2];

            return buffer;
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> point<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            if (static_cast<geometry_t>(serialized[0]) != geometry_t::POINT)
            {
                throw std::runtime_error(__tpf_error_message("Illegal identifier for deserialization of a point."));
            }

            return std::dynamic_pointer_cast<geometric_object<floatp_t>>(std::make_shared<point<floatp_t, kernel_t>>(
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 0 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 1 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 2 * sizeof(floatp_t)])
                ));
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Point_3& point<floatp_t, kernel_t>::get_internal() const
        {
            return this->_point;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::operator const typename kernel_t::Point_3&() const
        {
            return this->_point;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t>::operator typename kernel_t::Point_2() const
        {
            return typename kernel_t::Point_2(this->_point[0], this->_point[1]);
        }
    }
}
