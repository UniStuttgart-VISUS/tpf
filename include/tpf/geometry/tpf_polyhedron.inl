#include "tpf_polyhedron.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_tetrahedron.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../log/tpf_log.h"

#include "../stdext/tpf_comparator.h"

#include <CGAL/Delaunay_triangulation_3.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>::polyhedron(const std::vector<point<floatp_t, kernel_t>>& points) noexcept(false)
        {
            if (points.size() < 4)
            {
                throw std::runtime_error(__tpf_error_message("Number of points for a polyhedron must be at least 4."));
            }

            std::vector<typename CGAL::Delaunay_triangulation_3<kernel_t>::Point> delaunay_points;
            delaunay_points.reserve(points.size());

            std::transform(points.begin(), points.end(), std::back_inserter(delaunay_points), [](const point<floatp_t, kernel_t>& point)
            { return typename CGAL::Delaunay_triangulation_3<kernel_t>::Point(point.get_internal()); });

            CGAL::Delaunay_triangulation_3<kernel_t> delaunay(delaunay_points.begin(), delaunay_points.end());

            for (auto it = delaunay.finite_cells_begin(); it != delaunay.finite_cells_end(); ++it)
            {
                this->_tetrahedra.push_back(tetrahedron<floatp_t, kernel_t>(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point(), it->vertex(3)->point()));
            }

            if (this->_tetrahedra.size() == 0)
            {
                throw std::runtime_error(__tpf_error_message("Unable to create polyhedron. All points seem to be co-planar."));
            }
        }

        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>::polyhedron(const std::vector<tetrahedron<floatp_t, kernel_t>>& tetrahedra) noexcept
        {
            this->_tetrahedra = tetrahedra;
        }

        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>::polyhedron(const polyhedron<floatp_t, kernel_t>& copy) noexcept
        {
            this->_tetrahedra = copy._tetrahedra;
        }

        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>& polyhedron<floatp_t, kernel_t>::operator=(const std::vector<tetrahedron<floatp_t, kernel_t>>& tetrahedra) noexcept
        {
            geometric_object<floatp_t>::operator=(tetrahedra);

            this->_tetrahedra = tetrahedra;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>& polyhedron<floatp_t, kernel_t>::operator=(const polyhedron<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_tetrahedra = copy._tetrahedra;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool polyhedron<floatp_t, kernel_t>::operator==(const polyhedron<floatp_t, kernel_t>& other) const noexcept
        {
            return *this == other._tetrahedra;
        }

        template <typename floatp_t, typename kernel_t>
        inline bool polyhedron<floatp_t, kernel_t>::operator==(const std::vector<tetrahedron<floatp_t, kernel_t>>& other) const noexcept
        {
            // Pairwise compare points
            bool equal = this->_tetrahedra.size() == other.size();

            for (std::size_t o = 0; equal && o < this->_tetrahedra.size(); ++o)
            {
                for (const auto& outer_point : this->_tetrahedra[o].get_points())
                {
                    bool hit = false;

                    for (std::size_t i = 0; !hit && i < other.size(); ++i)
                    {
                        for (const auto& inner_point : other[i].get_points())
                        {
                            hit |= outer_point == inner_point;
                        }
                    }

                    equal &= hit;
                }
            }

            return equal;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> polyhedron<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<polyhedron<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& polyhedron<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                for (auto& tetrahedron : this->_tetrahedra)
                {
                    tetrahedron.transform(trafo);
                }
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_float<floatp_t, typename kernel_t::FT> polyhedron<floatp_t, kernel_t>::calculate_volume() const
        {
            auto volume = this->_tetrahedra[0].calculate_volume().get_exact_value();

            for (std::size_t i = 1; i < this->_tetrahedra.size(); ++i)
            {
                volume += this->_tetrahedra[i].calculate_volume().get_exact_value();
            }

            return volume;
        }

        template <typename floatp_t, typename kernel_t>
        inline point<floatp_t, kernel_t> polyhedron<floatp_t, kernel_t>::calculate_centroid() const
        {
            typename kernel_t::Vector_3 weighted_point(CGAL::NULL_VECTOR);
            typename kernel_t::FT weights = 0.0;

            for (const auto& tetrahedron : this->_tetrahedra)
            {
                const auto centroid = tetrahedron.calculate_centroid();
                const auto volume = tetrahedron.calculate_volume();

                weighted_point = weighted_point + volume.get_exact_value() * (centroid.get_internal() - CGAL::ORIGIN);
                weights += volume.get_exact_value();
            }

            weighted_point = weighted_point / weights;

            return CGAL::ORIGIN + weighted_point;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t polyhedron<floatp_t, kernel_t>::get_num_points() const
        {
            return this->get_points().size();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> polyhedron<floatp_t, kernel_t>::get_points() const
        {
            // Create set to get unique points
            std::set<Eigen::Matrix<floatp_t, 3, 1>> points_set;

            for (const auto& tetrahedron : this->_tetrahedra)
            {
                for (const auto& tetrahedron_point : tetrahedron.get_points())
                {
                    points_set.insert(tetrahedron_point);
                }
            }

            return std::vector<Eigen::Matrix<floatp_t, 3, 1>>(points_set.begin(), points_set.end());
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t polyhedron<floatp_t, kernel_t>::get_num_cells() const
        {
            return this->get_cells().size();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> polyhedron<floatp_t, kernel_t>::get_cells() const
        {
            // In get_points() we filter for unique points, therefore we need to map the cell indices to the correct points.
            // Build map of points to indices
            std::map<Eigen::Matrix<floatp_t, 3, 1>, std::size_t> point_idx_map;
            const auto &points = this->get_points();

            for (std::size_t i = 0; i < points.size(); ++i)
            {
                point_idx_map[points[i]] = i;
            }

            // Get all unique tetrahedron cells
            auto comparator = [](const std::vector<std::size_t>& lhs, const std::vector<std::size_t>& rhs)
            {
                return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
            };

            std::set<std::vector<std::size_t>, decltype(comparator)> cells(comparator);

            for (const auto& tetrahedron : this->_tetrahedra)
            {
                const auto tetrahedron_points = tetrahedron.get_points();
                auto tetrahedron_cells = tetrahedron.get_cells();

                for (auto& tetrahedron_cell : tetrahedron_cells)
                {
                    for (auto &tetrahedron_point_idx : tetrahedron_cell)
                    {
                        tetrahedron_point_idx = point_idx_map.at(tetrahedron_points[tetrahedron_point_idx]);
                    }

                    std::sort(tetrahedron_cell.begin(), tetrahedron_cell.end());

                    if (cells.find(tetrahedron_cell) != cells.end())
                    {
                        cells.erase(cells.find(tetrahedron_cell));
                    }
                    else
                    {
                        cells.insert(tetrahedron_cell);
                    }
                }
            }

            return std::vector<std::vector<std::size_t>>(cells.begin(), cells.end());
        }

        template <typename floatp_t, typename kernel_t>
        inline Eigen::Matrix<floatp_t, 3, 1> polyhedron<floatp_t, kernel_t>::get_centroid() const
        {
            return calculate_centroid().get_vertex();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> polyhedron<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> polyhedron<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline geometry_t polyhedron<floatp_t, kernel_t>::get_type() const
        {
            return geometry_t::POLYHEDRON;
        }

        template <typename floatp_t, typename kernel_t>
        inline const std::vector<tetrahedron<floatp_t, kernel_t>>& polyhedron<floatp_t, kernel_t>::get_internal() const
        {
            return this->_tetrahedra;
        }

        template <typename floatp_t, typename kernel_t>
        inline polyhedron<floatp_t, kernel_t>::operator const std::vector<tetrahedron<floatp_t, kernel_t>>&() const
        {
            return this->_tetrahedra;
        }
    }
}
