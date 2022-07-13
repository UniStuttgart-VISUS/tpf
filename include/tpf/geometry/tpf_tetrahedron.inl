#include "tpf_tetrahedron.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../log/tpf_log.h"

#include <CGAL/centroid.h>
#include <CGAL/number_utils.h>
#include <CGAL/Tetrahedron_3.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>::tetrahedron(const std::vector<point<floatp_t, kernel_t>>& points) noexcept(false)
        {
            if (points.size() != 4)
            {
                throw std::runtime_error(__tpf_error_message("Number of points for a tetrahedron must be 4."));
            }

            if (CGAL::coplanar<kernel_t>(points[0], points[1], points[2], points[3]))
            {
                throw std::runtime_error(__tpf_error_message("Points of a tetrahedron must not be coplanar."));
            }

            this->_tetrahedron = typename kernel_t::Tetrahedron_3(points[0], points[1], points[2], points[3]);
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>::tetrahedron(const point<floatp_t, kernel_t>& point_1,
            const point<floatp_t, kernel_t>& point_2, const point<floatp_t, kernel_t>& point_3, const point<floatp_t, kernel_t>& point_4) noexcept(false)
        {
            if (CGAL::coplanar<kernel_t>(point_1, point_2, point_3, point_4))
            {
                throw std::runtime_error(__tpf_error_message("Points of a tetrahedron must not be coplanar."));
            }

            this->_tetrahedron = typename kernel_t::Tetrahedron_3(point_1.get_internal(), point_2.get_internal(), point_3.get_internal(), point_4.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>::tetrahedron(const typename kernel_t::Tetrahedron_3& tetrahedron) noexcept
        {
            this->_tetrahedron = tetrahedron;
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>::tetrahedron(const tetrahedron<floatp_t, kernel_t>& copy) noexcept
        {
            this->_tetrahedron = copy._tetrahedron;
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>& tetrahedron<floatp_t, kernel_t>::operator=(const typename kernel_t::Tetrahedron_3& tetrahedron) noexcept
        {
            geometric_object<floatp_t>::operator=(tetrahedron);

            this->_tetrahedron = tetrahedron;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>& tetrahedron<floatp_t, kernel_t>::operator=(const tetrahedron<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_tetrahedron = copy._tetrahedron;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool tetrahedron<floatp_t, kernel_t>::operator==(const tetrahedron<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_tetrahedron == other._tetrahedron;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool tetrahedron<floatp_t, kernel_t>::operator==(const typename kernel_t::Tetrahedron_3& other) const noexcept
        {
            return this->_tetrahedron == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> tetrahedron<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<tetrahedron<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& tetrahedron<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                this->_tetrahedron = typename kernel_t::Tetrahedron_3(
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(0)).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(1)).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(2)).transform(trafo)).get_internal(),
                    static_cast<point<floatp_t, kernel_t>&>(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(3)).transform(trafo)).get_internal());
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> tetrahedron<floatp_t, kernel_t>::calculate_volume() const
        {
            return CGAL::abs(this->_tetrahedron.volume());
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> tetrahedron<floatp_t, kernel_t>::calculate_centroid() const
        {
            return CGAL::centroid(this->_tetrahedron);
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t tetrahedron<floatp_t, kernel_t>::get_num_points() const
        {
            return 4;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> tetrahedron<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.push_back(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(0)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(1)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(2)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_tetrahedron.vertex(3)).get_vertex());

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t tetrahedron<floatp_t, kernel_t>::get_num_cells() const
        {
            return 4;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> tetrahedron<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(4));

            // TODO winding order of vertices?
            cells[0].push_back(static_cast<std::size_t>(2));
            cells[0].push_back(static_cast<std::size_t>(1));
            cells[0].push_back(static_cast<std::size_t>(0));

            cells[1].push_back(static_cast<std::size_t>(3));
            cells[1].push_back(static_cast<std::size_t>(1));
            cells[1].push_back(static_cast<std::size_t>(0));

            cells[2].push_back(static_cast<std::size_t>(3));
            cells[2].push_back(static_cast<std::size_t>(2));
            cells[2].push_back(static_cast<std::size_t>(0));

            cells[3].push_back(static_cast<std::size_t>(3));
            cells[3].push_back(static_cast<std::size_t>(2));
            cells[3].push_back(static_cast<std::size_t>(1));

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> tetrahedron<floatp_t, kernel_t>::get_centroid() const
        {
            return (1.0 / 4.0) * (point<floatp_t, kernel_t>(this->_tetrahedron.vertex(0)).get_vertex()
                + point<floatp_t, kernel_t>(this->_tetrahedron.vertex(1)).get_vertex()
                + point<floatp_t, kernel_t>(this->_tetrahedron.vertex(2)).get_vertex()
                + point<floatp_t, kernel_t>(this->_tetrahedron.vertex(3)).get_vertex());
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> tetrahedron<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> tetrahedron<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline geometry_t tetrahedron<floatp_t, kernel_t>::get_type() const
        {
            return geometry_t::TETRAHEDRON;
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Tetrahedron_3& tetrahedron<floatp_t, kernel_t>::get_internal() const
        {
            return this->_tetrahedron;
        }

        template <typename floatp_t, typename kernel_t>
        inline tetrahedron<floatp_t, kernel_t>::operator const typename kernel_t::Tetrahedron_3&() const
        {
            return this->_tetrahedron;
        }
    }
}
