#include "tpf_cuboid.h"

#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../exception/tpf_not_implemented_exception.h"

#include <CGAL/Iso_cuboid_3.h>

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>::cuboid(const point<floatp_t, kernel_t>& point_min, const point<floatp_t, kernel_t>& point_max) noexcept
        {
            this->_cuboid = typename kernel_t::Iso_cuboid_3(point_min.get_internal(), point_max.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>::cuboid(const typename kernel_t::Iso_cuboid_3& cuboid) noexcept
        {
            this->_cuboid = cuboid;
        }

        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>::cuboid(const cuboid<floatp_t, kernel_t>& copy) noexcept
        {
            this->_cuboid = copy._cuboid;
        }

        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>& cuboid<floatp_t, kernel_t>::operator=(const typename kernel_t::Iso_cuboid_3& cuboid) noexcept
        {
            geometric_object<floatp_t>::operator=(cuboid);

            this->_cuboid = cuboid;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>& cuboid<floatp_t, kernel_t>::operator=(const cuboid<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_cuboid = copy._cuboid;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool cuboid<floatp_t, kernel_t>::operator==(const cuboid<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_cuboid == other._cuboid;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool cuboid<floatp_t, kernel_t>::operator==(const typename kernel_t::Iso_cuboid_3& other) const noexcept
        {
            return this->_cuboid == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> cuboid<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<cuboid<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& cuboid<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                this->_cuboid = typename kernel_t::Iso_cuboid_3(
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_cuboid.min()).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_cuboid.max()).transform(trafo)).get_internal());
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> cuboid<floatp_t, kernel_t>::calculate_volume() const
        {
            return this->_cuboid.volume();
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> cuboid<floatp_t, kernel_t>::get_min_point() const
        {
            return static_cast<point<floatp_t, kernel_t>>(this->_cuboid.min());
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> cuboid<floatp_t, kernel_t>::get_max_point() const
        {
            return static_cast<point<floatp_t, kernel_t>>(this->_cuboid.max());
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> cuboid<floatp_t, kernel_t>::get_center_point() const
        {
            return this->_cuboid.min() + 0.5 * (this->_cuboid.max() - this->_cuboid.min());
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t cuboid<floatp_t, kernel_t>::get_num_points() const
        {
            return 8;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> cuboid<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(0)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(1)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(2)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(3)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(4)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(5)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(6)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_cuboid.vertex(7)).get_vertex());

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t cuboid<floatp_t, kernel_t>::get_num_cells() const
        {
            return 6;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> cuboid<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(6));

            // vertex id to coords mapping:
            // vert. | x | y | z
            // ------+---+---+---
            //     0 | 0 | 0 | 0
            //     1 | 1 | 0 | 0
            //     2 | 1 | 1 | 0
            //     3 | 0 | 1 | 0
            //     4 | 0 | 1 | 1
            //     5 | 0 | 0 | 1
            //     6 | 1 | 0 | 1
            //     7 | 1 | 1 | 1

            // draw all quads counter-clockwise when viewed from outside

            // z = 0 plane
            cells[0].push_back(static_cast<std::size_t>(0));
            cells[0].push_back(static_cast<std::size_t>(3));
            cells[0].push_back(static_cast<std::size_t>(2));
            cells[0].push_back(static_cast<std::size_t>(1));

            // z = 1 plane
            cells[1].push_back(static_cast<std::size_t>(5));
            cells[1].push_back(static_cast<std::size_t>(6));
            cells[1].push_back(static_cast<std::size_t>(7));
            cells[1].push_back(static_cast<std::size_t>(4));

            // y = 0 plane
            cells[2].push_back(static_cast<std::size_t>(0));
            cells[2].push_back(static_cast<std::size_t>(1));
            cells[2].push_back(static_cast<std::size_t>(6));
            cells[2].push_back(static_cast<std::size_t>(5));

            // y = 1 plane
            cells[3].push_back(static_cast<std::size_t>(3));
            cells[3].push_back(static_cast<std::size_t>(4));
            cells[3].push_back(static_cast<std::size_t>(7));
            cells[3].push_back(static_cast<std::size_t>(2));

            // x = 0 plane
            cells[4].push_back(static_cast<std::size_t>(0));
            cells[4].push_back(static_cast<std::size_t>(5));
            cells[4].push_back(static_cast<std::size_t>(4));
            cells[4].push_back(static_cast<std::size_t>(3));

            // x = 1 plane
            cells[5].push_back(static_cast<std::size_t>(1));
            cells[5].push_back(static_cast<std::size_t>(2));
            cells[5].push_back(static_cast<std::size_t>(7));
            cells[5].push_back(static_cast<std::size_t>(6));

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> cuboid<floatp_t, kernel_t>::get_centroid() const
        {
            return (1.0 / 8.0) * (point<floatp_t, kernel_t>(this->_cuboid.vertex(0)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(1)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(2)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(3)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(4)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(5)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(6)).get_vertex()
                + point<floatp_t, kernel_t>(this->_cuboid.vertex(7)).get_vertex());
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> cuboid<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> cuboid<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Iso_cuboid_3& cuboid<floatp_t, kernel_t>::get_internal() const
        {
            return this->_cuboid;
        }

        template <typename floatp_t, typename kernel_t>
        inline cuboid<floatp_t, kernel_t>::operator const typename kernel_t::Iso_cuboid_3&() const
        {
            return this->_cuboid;
        }
    }
}
