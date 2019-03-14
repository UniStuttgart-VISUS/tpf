#include "tpf_intersection.h"

#include "tpf_cuboid.h"
#include "tpf_line.h"
#include "tpf_plane.h"
#include "tpf_point.h"
#include "tpf_polyhedron.h"
#include "tpf_tetrahedron.h"
#include "tpf_triangle.h"

#include "../algorithm/tpf_jooat.h"

#include "../stdext/tpf_comparator.h"

#include "../utility/tpf_optional.h"

#include "Eigen/Dense"

#include <boost/variant/get.hpp>

#include <CGAL/intersections.h>
#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline bool does_intersect_with(const line<floatp_t, kernel_t>& line, const plane<floatp_t, kernel_t>& plane)
        {
            return CGAL::do_intersect(plane.get_internal(), line.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline bool does_intersect_with(const plane<floatp_t, kernel_t>& plane, const cuboid<floatp_t, kernel_t>& cuboid)
        {
            return CGAL::do_intersect(cuboid.get_internal().bbox(), plane.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline utility::optional<point<floatp_t, kernel_t>> intersect_with(const line<floatp_t, kernel_t>& line, const plane<floatp_t, kernel_t>& plane)
        {
            auto intersection = CGAL::intersection(plane.get_internal(), line.get_internal());

            if (intersection)
            {
                const typename kernel_t::Point_3* p = boost::get<typename kernel_t::Point_3>(&*intersection);

                if (p != nullptr)
                {
                    return *p;
                }
                else
                {
                    return utility::nullopt;
                }
            }
            else
            {
                return utility::nullopt;
            }
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const cuboid<floatp_t, kernel_t>& cuboid)
        {
            // Extract edges
            std::vector<line<floatp_t, kernel_t>> edges;
            edges.reserve(12);

            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(0), cuboid.get_internal().vertex(1)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(3), cuboid.get_internal().vertex(2)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(5), cuboid.get_internal().vertex(6)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(4), cuboid.get_internal().vertex(7)));

            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(0), cuboid.get_internal().vertex(3)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(1), cuboid.get_internal().vertex(2)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(5), cuboid.get_internal().vertex(4)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(6), cuboid.get_internal().vertex(7)));

            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(0), cuboid.get_internal().vertex(5)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(1), cuboid.get_internal().vertex(6)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(3), cuboid.get_internal().vertex(4)));
            edges.push_back(line<floatp_t, kernel_t>(cuboid.get_internal().vertex(2), cuboid.get_internal().vertex(7)));

            // Intersect edges with plane
            std::vector<point<floatp_t, kernel_t>> intersections;

            for (const line<floatp_t, kernel_t>& edge : edges)
            {
                if (does_intersect_with(edge, plane))
                {
                    intersections.push_back(*intersect_with(edge, plane));
                }
            }

            return intersections;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const tetrahedron<floatp_t, kernel_t>& tetrahedron)
        {
            // Extract edges
            std::vector<line<floatp_t, kernel_t>> edges;
            edges.reserve(6);

            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(0), tetrahedron.get_internal().vertex(1)));
            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(0), tetrahedron.get_internal().vertex(2)));
            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(0), tetrahedron.get_internal().vertex(3)));
            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(1), tetrahedron.get_internal().vertex(2)));
            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(1), tetrahedron.get_internal().vertex(3)));
            edges.push_back(line<floatp_t, kernel_t>(tetrahedron.get_internal().vertex(2), tetrahedron.get_internal().vertex(3)));
            
            // Intersect edges with plane
            std::vector<point<floatp_t, kernel_t>> intersections;

            for (const line<floatp_t, kernel_t>& edge : edges)
            {
                if (does_intersect_with(edge, plane))
                {
                    intersections.push_back(*intersect_with(edge, plane));
                }
            }

            return intersections;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<point<floatp_t, kernel_t>> intersect_with(const plane<floatp_t, kernel_t>& plane, const polyhedron<floatp_t, kernel_t>& polyhedron)
        {
            // Extract faces
            auto face_predicate = [](const triangle<floatp_t, kernel_t>& triangle) -> std::size_t
            {
                auto points = triangle.get_points();
                std::sort(points.begin(), points.end(), std::less<Eigen::Matrix<floatp_t, 3, 1>>());

                return algorithm::joaat_hash(points[0], points[1], points[2]);
            };

            std::unordered_map<triangle<floatp_t, kernel_t>, std::size_t, decltype(face_predicate)> faces(23, face_predicate);

            for (const auto& tetrahedron : polyhedron.get_internal())
            {
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = i + 1; j < 3; ++j)
                    {
                        for (int k = j + 1; k < 4; ++k)
                        {
                            const triangle<floatp_t, kernel_t> face(tetrahedron.get_internal().vertex(i),
                                tetrahedron.get_internal().vertex(j), tetrahedron.get_internal().vertex(k));

                            if (faces.find(face) == faces.end())
                            {
                                faces[face] = 1;
                            }
                            else
                            {
                                ++faces[face];
                            }
                        }
                    }
                }
            }

            // Filter faces, such that only outer ones remain
            for (auto it = faces.begin(); it != faces.end(); )
            {
                if (it->second != 1)
                {
                    faces.erase(it++);
                }
                else
                {
                    ++it;
                }
            }

            // Extract edges
            auto edge_predicate = [](const line<floatp_t, kernel_t>& line) -> std::size_t
            {
                auto points = line.get_points();
                std::sort(points.begin(), points.end(), std::less<Eigen::Matrix<floatp_t, 3, 1>>());

                return algorithm::joaat_hash(points[0], points[1]);
            };

            std::unordered_set<line<floatp_t, kernel_t>, decltype(edge_predicate)> edges(23, edge_predicate);

            for (const auto& face : faces)
            {
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = i + 1; j < 3; ++j)
                    {
                        const line<floatp_t, kernel_t> edge(face.first.get_internal().vertex(i), face.first.get_internal().vertex(j));

                        if (edges.find(edge) == edges.end())
                        {
                            edges.insert(edge);
                        }
                    }
                }
            }

            // Intersect edges with plane
            std::vector<point<floatp_t, kernel_t>> intersections;

            for (const auto& edge : edges)
            {
                if (does_intersect_with(edge, plane))
                {
                    intersections.push_back(*intersect_with(edge, plane));
                }
            }

            return intersections;
        }
    }
}