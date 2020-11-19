#include "tpf_plane.h"

#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../math/tpf_vector.h"

#include "Eigen/Dense"

#include <CGAL/Plane_3.h>

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>::plane(const point<floatp_t, kernel_t>& first, const Eigen::Matrix<floatp_t, 3, 1>& normal) noexcept
        {
            const auto vectors = math::template orthonormal<floatp_t>(normal);

            const point<floatp_t, kernel_t> second(first.get_vertex() + vectors.first);
            const point<floatp_t, kernel_t> third(first.get_vertex() + vectors.second);

            this->_plane = typename kernel_t::Plane_3(first.get_internal(), second.get_internal(), third.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>::plane(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2, const point<floatp_t, kernel_t>& point_3) noexcept
        {
            this->_plane = typename kernel_t::Plane_3(point_1.get_internal(), point_2.get_internal(), point_3.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>::plane(const typename kernel_t::Plane_3& plane) noexcept
        {
            this->_plane = plane;
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>::plane(const plane<floatp_t, kernel_t>& copy) noexcept
        {
            this->_plane = copy._plane;
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>& plane<floatp_t, kernel_t>::operator=(const typename kernel_t::Plane_3& plane) noexcept
        {
            geometric_object<floatp_t>::operator=(plane);

            this->_plane = plane;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>& plane<floatp_t, kernel_t>::operator=(const plane<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_plane = copy._plane;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool plane<floatp_t, kernel_t>::operator==(const plane<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_plane == other._plane;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool plane<floatp_t, kernel_t>::operator==(const typename kernel_t::Plane_3& other) const noexcept
        {
            return this->_plane == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> plane<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<plane<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& plane<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                this->_plane = typename kernel_t::Plane_3(
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_plane.point()).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_plane.point() + this->_plane.base1()).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_plane.point() + this->_plane.base2()).transform(trafo)).get_internal());
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t plane<floatp_t, kernel_t>::get_num_points() const
        {
            return 0;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> plane<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t plane<floatp_t, kernel_t>::get_num_cells() const
        {
            return 0;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> plane<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells;
            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> plane<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> plane<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Plane_3& plane<floatp_t, kernel_t>::get_internal() const
        {
            return this->_plane;
        }

        template <typename floatp_t, typename kernel_t>
        inline plane<floatp_t, kernel_t>::operator const typename kernel_t::Plane_3&() const
        {
            return this->_plane;
        }
    }
}
