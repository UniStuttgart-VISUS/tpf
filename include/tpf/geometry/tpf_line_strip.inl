#include "tpf_line_strip.h"

#include "tpf_geometric_object.h"
#include "tpf_line.h"
#include "tpf_point.h"
#include "tpf_polyhedron.h"
#include "tpf_triangle.h"

#include "../log/tpf_log.h"

#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline line_strip<floatp_t, kernel_t>::line_strip(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2) noexcept
        {
            this->_line_strip.push_back(point_1.get_internal());
            this->_line_strip.push_back(point_2.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline line_strip<floatp_t, kernel_t>::line_strip(const std::vector<point<floatp_t, kernel_t>>& points) noexcept
        {
            this->_line_strip = points;
        }

        template <typename floatp_t, typename kernel_t>
        inline line_strip<floatp_t, kernel_t>::line_strip(const line_strip<floatp_t, kernel_t>& copy) noexcept
        {
            this->_line_strip = copy._line_strip;
        }

        template <typename floatp_t, typename kernel_t>
        inline line_strip<floatp_t, kernel_t>& line_strip<floatp_t, kernel_t>::operator=(const line_strip<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_line_strip = copy._line_strip;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool line_strip<floatp_t, kernel_t>::operator==(const line_strip<floatp_t, kernel_t>& other) const noexcept
        {
            bool equal = this->_line_strip.size() == other._line_strip.size();

            for (std::size_t i = 0; i < this->_line_strip.size() && equal; ++i)
            {
                equal &= (this->_line_strip[i] == other._line_strip[i]);
            }

            return equal;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> line_strip<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<line_strip<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& line_strip<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                for (auto& point : this->_line_strip)
                {
                    point.transform(trafo);
                }
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t line_strip<floatp_t, kernel_t>::get_num_points() const
        {
            return this->_line_strip.size();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> line_strip<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.reserve(this->_line_strip.size());

            for (const auto& point : this->_line_strip)
            {
                points.push_back(point.get_vertex());
            }

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t line_strip<floatp_t, kernel_t>::get_num_cells() const
        {
            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> line_strip<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(1));
            cells.front().resize(this->_line_strip.size());

            std::iota(cells.front().begin(), cells.front().end(), 0);

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> line_strip<floatp_t, kernel_t>::get_centroid() const
        {
            switch (this->_line_strip.size())
            {
            case 2:
                return line<floatp_t, kernel_t>(this->_line_strip[0], this->_line_strip[1]).get_centroid();
            case 3:
                return triangle<floatp_t, kernel_t>(this->_line_strip[0], this->_line_strip[1], this->_line_strip[2]).get_centroid();
            default:
                return polyhedron<floatp_t, kernel_t>(this->_line_strip).calculate_centroid();
            }
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> line_strip<floatp_t, kernel_t>::serialize() const
        {
            std::vector<char> buffer(sizeof(uint64_t) + sizeof(uint64_t) + 3 * this->_line_strip.size() * sizeof(floatp_t));

            *reinterpret_cast<uint64_t*>(&buffer[0]) = static_cast<uint64_t>(get_type());

            *reinterpret_cast<uint64_t*>(&buffer[sizeof(uint64_t)]) = static_cast<uint64_t>(this->_line_strip.size());

            for (std::size_t i = 0; i < this->_line_strip.size(); ++i)
            {
                const auto offset = sizeof(uint64_t) + sizeof(uint64_t) + 3 * i * sizeof(floatp_t);

                *reinterpret_cast<floatp_t*>(&buffer[offset + 0 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line_strip[i].get_internal().x()));
                *reinterpret_cast<floatp_t*>(&buffer[offset + 1 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line_strip[i].get_internal().y()));
                *reinterpret_cast<floatp_t*>(&buffer[offset + 2 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line_strip[i].get_internal().z()));
            }

            return buffer;
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> line_strip<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            if (static_cast<geometry_t>(*reinterpret_cast<const uint64_t*>(&serialized[0])) != get_type())
            {
                throw std::runtime_error(__tpf_error_message("Illegal identifier for deserialization of a line strip."));
            }

            const auto num_points = *reinterpret_cast<const uint64_t*>(&serialized[sizeof(uint64_t)])

            std::vector<point<floatp_t, kernel_t>> points;
            points.reserve(num_points);

            auto offset = sizeof(uint64_t) + sizeof(uint64_t);

            for (std::size_t i = 0; i < num_points; ++i)
            {
                points.push_back(point<floatp_t, kernel_t>(
                    *reinterpret_cast<const floatp_t*>(&serialized[offset + 0 * sizeof(floatp_t)]),
                    *reinterpret_cast<const floatp_t*>(&serialized[offset + 1 * sizeof(floatp_t)]),
                    *reinterpret_cast<const floatp_t*>(&serialized[offset + 2 * sizeof(floatp_t)])));

                offset += 3 * sizeof(floatp_t);
            }

            return std::dynamic_pointer_cast<geometric_object<floatp_t>>(std::make_shared<line_strip<floatp_t, kernel_t>>(points));
        }

        template <typename floatp_t, typename kernel_t>
        inline geometry_t line_strip<floatp_t, kernel_t>::get_type() const
        {
            return geometry_t::LINE_STRIP;
        }

        template <typename floatp_t, typename kernel_t>
        inline const std::vector<point<floatp_t, kernel_t>>& line_strip<floatp_t, kernel_t>::get_internal() const
        {
            return this->_line_strip;
        }

        template <typename floatp_t, typename kernel_t>
        inline line_strip<floatp_t, kernel_t>::operator const std::vector<point<floatp_t, kernel_t>>&() const
        {
            return this->_line_strip;
        }
    }
}
