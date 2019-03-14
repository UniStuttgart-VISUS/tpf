#include "tpf_line.h"

#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "../log/tpf_log.h"

#include <CGAL/Line_3.h>

#include <memory>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>::line(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2) noexcept
        {
            this->_line = typename kernel_t::Segment_3(point_1.get_internal(), point_2.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>::line(const typename kernel_t::Segment_3& line) noexcept
        {
            this->_line = line;
        }

        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>::line(const line<floatp_t, kernel_t>& copy) noexcept
        {
            this->_line = copy._line;
        }

        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>& line<floatp_t, kernel_t>::operator=(const typename kernel_t::Segment_3& line) noexcept
        {
            geometric_object<floatp_t>::operator=(line);

            this->_line = line;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>& line<floatp_t, kernel_t>::operator=(const line<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_line = copy._line;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool line<floatp_t, kernel_t>::operator==(const line<floatp_t, kernel_t>& other) const noexcept
        {
            return this->_line == other._line;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool line<floatp_t, kernel_t>::operator==(const typename kernel_t::Segment_3& other) const noexcept
        {
            return this->_line == other;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> line<floatp_t, kernel_t>::clone() const
        {
            return std::make_shared<line<floatp_t, kernel_t>>(*this);
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t line<floatp_t, kernel_t>::get_num_points() const
        {
            return 2;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> line<floatp_t, kernel_t>::get_points() const
        {
            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.push_back(point<floatp_t, kernel_t>(this->_line.vertex(0)).get_vertex());
            points.push_back(point<floatp_t, kernel_t>(this->_line.vertex(1)).get_vertex());
            
            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t line<floatp_t, kernel_t>::get_num_cells() const
        {
            return 1;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> line<floatp_t, kernel_t>::get_cells() const
        {
            std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(1));
            cells.front().push_back(static_cast<std::size_t>(0));
            cells.front().push_back(static_cast<std::size_t>(1));

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> line<floatp_t, kernel_t>::serialize() const
        {
            std::vector<char> buffer(1 + 6 * sizeof(floatp_t));

            buffer[0] = static_cast<char>(geometry_t::LINE);

            *reinterpret_cast<floatp_t*>(&buffer[1 + 0 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(0).x()));
            *reinterpret_cast<floatp_t*>(&buffer[1 + 1 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(0).y()));
            *reinterpret_cast<floatp_t*>(&buffer[1 + 2 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(0).z()));
            *reinterpret_cast<floatp_t*>(&buffer[1 + 3 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(1).x()));
            *reinterpret_cast<floatp_t*>(&buffer[1 + 4 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(1).y()));
            *reinterpret_cast<floatp_t*>(&buffer[1 + 5 * sizeof(floatp_t)]) = static_cast<floatp_t>(CGAL::to_double(this->_line.vertex(1).z()));

            return buffer;
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> line<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            if (static_cast<geometry_t>(serialized[0]) != geometry_t::LINE)
            {
                throw std::runtime_error(__tpf_error_message("Illegal identifier for deserialization of a line."));
            }

            return std::dynamic_pointer_cast<geometric_object<floatp_t>>(std::make_shared<line<floatp_t, kernel_t>>(point<floatp_t, kernel_t>(
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 0 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 1 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 2 * sizeof(floatp_t)])), point<floatp_t, kernel_t>(
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 3 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 4 * sizeof(floatp_t)]),
                *reinterpret_cast<const floatp_t*>(&serialized[1 + 5 * sizeof(floatp_t)]))
                ));
        }

        template <typename floatp_t, typename kernel_t>
        inline const typename kernel_t::Segment_3& line<floatp_t, kernel_t>::get_internal() const
        {
            return this->_line;
        }

        template <typename floatp_t, typename kernel_t>
        inline line<floatp_t, kernel_t>::operator const typename kernel_t::Segment_3&() const
        {
            return this->_line;
        }
    }
}
