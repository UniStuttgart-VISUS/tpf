#include "tpf_factory.h"

#include "tpf_geometric_object.h"

#include "tpf_cuboid.h"
#include "tpf_line.h"
#include "tpf_plane.h"
#include "tpf_point.h"
#include "tpf_polygon.h"
#include "tpf_polyhedron.h"
#include "tpf_rectangle.h"
#include "tpf_tetrahedron.h"
#include "tpf_triangle.h"

#include "../log/tpf_log.h"

#include "../utility/tpf_optional.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kernel/global_functions_3.h>

#include "Eigen/Dense"

#include <algorithm>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        namespace
        {
            /// <summary>
            /// Default make routine for geometric objects that can be safely created by calling the constructor directly
            /// </summary>
            /// <template name="object_t">Geometric object type</template>
            /// <template name="non_static_make_shared">Indicate that this was called from the non-static make_shared function</template>
            template <typename object_t, bool non_static_make_shared>
            struct make_impl
            {
                /// <summary>
                /// Create geometric object by calling its constructor with the given arguments
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_nothrow_constructible<object_t, arguments_t...>::value, utility::optional<object_t>>::type make(arguments_t... args) noexcept
                {
                    return object_t(args...);
                }

                /// <summary>
                /// Return unfilled optional if the object is not constructible
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                template <typename... arguments_t>
                static typename std::enable_if<!std::is_constructible<object_t, arguments_t...>::value, utility::optional<object_t>>::type make(arguments_t...) noexcept
                {
                    return utility::nullopt;
                }

                /// <summary>
                /// Fail, if there is only a throwing constructor
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_constructible<object_t, arguments_t...>::value && !std::is_nothrow_constructible<object_t, arguments_t...>::value,
                    utility::optional<object_t>>::type make(arguments_t...) noexcept
                {
                    static_assert(non_static_make_shared, "For this constructor, a safe factory-make is missing.");
                    return utility::nullopt;
                }
            };

            /// <summary>
            /// Specialization for tetrahedra
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            /// <template name="kernel_t">CGAL kernel</template>
            template <typename floatp_t, typename kernel_t, bool non_static_make_shared>
            struct make_impl<tetrahedron<floatp_t, kernel_t>, non_static_make_shared>
            {
                /// <summary>
                /// Create geometric object using a safe constructor with the given arguments
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_nothrow_constructible<tetrahedron<floatp_t, kernel_t>, arguments_t...>::value,
                    utility::optional<tetrahedron<floatp_t, kernel_t>>>::type make(arguments_t... args) noexcept
                {
                    return tetrahedron<floatp_t, kernel_t>(args...);
                }

                /// <summary>
                /// Return unfilled optional if the object is not constructible
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<!std::is_constructible<tetrahedron<floatp_t, kernel_t>, arguments_t...>::value, utility::optional<tetrahedron<floatp_t, kernel_t>>>::type
                    make(arguments_t...) noexcept
                {
                    return utility::nullopt;
                }

                /// <summary>
                /// Fail, if there is only a throwing constructor
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_constructible<tetrahedron<floatp_t, kernel_t>, arguments_t...>::value &&
                    !std::is_nothrow_constructible<tetrahedron<floatp_t, kernel_t>, arguments_t...>::value,
                    utility::optional<tetrahedron<floatp_t, kernel_t>>>::type make(arguments_t...) noexcept
                {
                    static_assert(non_static_make_shared, "For this constructor, a safe factory-make is missing.");
                    return utility::nullopt;
                }

                /// <summary>
                /// Create geometric object by error-checking and then calling its constructor with the given arguments
                /// </summary>
                static utility::optional<tetrahedron<floatp_t, kernel_t>> make(const std::vector<point<floatp_t, kernel_t>>& points) noexcept
                {
                    return (points.size() != 4 || CGAL::template coplanar<kernel_t>(points[0], points[1], points[2], points[3]))
                        ? utility::optional<tetrahedron<floatp_t, kernel_t>>(utility::nullopt) : tetrahedron<floatp_t, kernel_t>(points);
                }

                /// <summary>
                /// Create geometric object by error-checking and then calling its constructor with the given arguments
                /// </summary>
                static utility::optional<tetrahedron<floatp_t, kernel_t>> make(const point<floatp_t, kernel_t>& point_1,
                    const point<floatp_t, kernel_t>& point_2, const point<floatp_t, kernel_t>& point_3, const point<floatp_t, kernel_t>& point_4) noexcept
                {
                    return CGAL::template coplanar<kernel_t>(point_1, point_2, point_3, point_4)
                        ? utility::optional<tetrahedron<floatp_t, kernel_t>>(utility::nullopt) : tetrahedron<floatp_t, kernel_t>(point_1, point_2, point_3, point_4);
                }
            };

            /// <summary>
            /// Specialization for polygons
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            /// <template name="kernel_t">CGAL kernel</template>
            template <typename floatp_t, typename kernel_t, bool non_static_make_shared>
            struct make_impl<polygon<floatp_t, kernel_t>, non_static_make_shared>
            {
                /// <summary>
                /// Create geometric object using a safe constructor with the given arguments
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_nothrow_constructible<polygon<floatp_t, kernel_t>, arguments_t...>::value, utility::optional<polygon<floatp_t, kernel_t>>>::type
                    make(arguments_t... args) noexcept
                {
                    return polygon<floatp_t, kernel_t>(args...);
                }

                /// <summary>
                /// Return unfilled optional if the object is not constructible
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<!std::is_constructible<polygon<floatp_t, kernel_t>, arguments_t...>::value, utility::optional<polygon<floatp_t, kernel_t>>>::type
                    make(arguments_t...) noexcept
                {
                    return utility::nullopt;
                }

                /// <summary>
                /// Fail, if there is only a throwing constructor
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_constructible<polygon<floatp_t, kernel_t>, arguments_t...>::value &&
                    !std::is_nothrow_constructible<polygon<floatp_t, kernel_t>, arguments_t...>::value,
                    utility::optional<polygon<floatp_t, kernel_t>>>::type make(arguments_t...) noexcept
                {
                    static_assert(non_static_make_shared, "For this constructor, a safe factory-make is missing.");
                    return utility::nullopt;
                }

                /// <summary>
                /// Create geometric object by error-checking and then calling its constructor with the given arguments
                /// </summary>
                static utility::optional<polygon<floatp_t, kernel_t>> make(const std::vector<point<floatp_t, kernel_t>>& points,
                    const Eigen::Matrix<floatp_t, 3, 1>& normal, const bool make_convex, const bool render_as_triangles = false) noexcept
                {
                    return (points.size() < 3) ? utility::optional<polygon<floatp_t, kernel_t>>(utility::nullopt)
                        : polygon<floatp_t, kernel_t>(points, normal, make_convex, render_as_triangles);
                }

                /// <summary>
                /// Create geometric object by error-checking and then calling its constructor with the given arguments
                /// </summary>
                static utility::optional<polygon<floatp_t, kernel_t>> make(const std::vector<point<floatp_t, kernel_t>>& points,
                    const bool make_convex, const bool render_as_triangles = false) noexcept
                {
                    return (points.size() < 3) ? utility::optional<polygon<floatp_t, kernel_t>>(utility::nullopt)
                        : polygon<floatp_t, kernel_t>(points, make_convex, render_as_triangles);
                }
            };

            /// <summary>
            /// Specialization for polyhedra
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            /// <template name="kernel_t">CGAL kernel</template>
            template <typename floatp_t, typename kernel_t, bool non_static_make_shared>
            struct make_impl<polyhedron<floatp_t, kernel_t>, non_static_make_shared>
            {
                /// <summary>
                /// Create geometric object using a safe constructor with the given arguments
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_nothrow_constructible<polyhedron<floatp_t, kernel_t>, arguments_t...>::value, utility::optional<polyhedron<floatp_t, kernel_t>>>::type
                    make(arguments_t... args) noexcept
                {
                    return polyhedron<floatp_t, kernel_t>(args...);
                }

                /// <summary>
                /// Return unfilled optional if the object is not constructible
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                /// <param name="args">Constructor parameters</param>
                template <typename... arguments_t>
                static typename std::enable_if<!std::is_constructible<polyhedron<floatp_t, kernel_t>, arguments_t...>::value, utility::optional<polyhedron<floatp_t, kernel_t>>>::type
                    make(arguments_t...) noexcept
                {
                    return utility::nullopt;
                }

                /// <summary>
                /// Fail, if there is only a throwing constructor
                /// </summary>
                /// <template name="arguments_t">Type of constructor parameters</template>
                template <typename... arguments_t>
                static typename std::enable_if<std::is_constructible<polyhedron<floatp_t, kernel_t>, arguments_t...>::value &&
                    !std::is_nothrow_constructible<polyhedron<floatp_t, kernel_t>, arguments_t...>::value,
                    utility::optional<polyhedron<floatp_t, kernel_t>>>::type make(arguments_t...) noexcept
                {
                    static_assert(non_static_make_shared, "For this constructor, a safe factory-make is missing.");
                    return utility::nullopt;
                }

                /// <summary>
                /// Create geometric object by error-checking and then calling its constructor with the given arguments
                /// </summary>
                static utility::optional<polyhedron<floatp_t, kernel_t>> make(const std::vector<point<floatp_t, kernel_t>>& points) noexcept
                {
                    if (points.size() < 4)
                    {
                        return utility::nullopt;
                    }

                    std::vector<tetrahedron<floatp_t, kernel_t>> tetrahedra;

                    std::vector<typename CGAL::Delaunay_triangulation_3<kernel_t>::Point> delaunay_points;
                    delaunay_points.reserve(points.size());

                    std::transform(points.begin(), points.end(), std::back_inserter(delaunay_points), [](const point<floatp_t, kernel_t>& point)
                    { return typename CGAL::Delaunay_triangulation_3<kernel_t>::Point(point.get_internal()); });

                    CGAL::Delaunay_triangulation_3<kernel_t> delaunay(delaunay_points.begin(), delaunay_points.end());

                    for (auto it = delaunay.finite_cells_begin(); it != delaunay.finite_cells_end(); ++it)
                    {
                        tetrahedra.push_back(tetrahedron<floatp_t, kernel_t>(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point(), it->vertex(3)->point()));
                    }

                    return (tetrahedra.size() == 0) ? utility::optional<polyhedron<floatp_t, kernel_t>>(utility::nullopt) : polyhedron<floatp_t, kernel_t>(tetrahedra);
                }
            };

            /// <summary>
            /// Default make routine for geometric objects that can be safely created by calling the constructor directly
            /// </summary>
            /// <template name="object_t">Geometric object type</template>
            template <typename object_t>
            inline std::shared_ptr<object_t> make_shared(utility::optional<object_t>&& object)
            {
                return object ? std::make_shared<object_t>(*object) : nullptr;
            }
        }

        template <typename object_t, typename... arguments_t>
        static inline utility::optional<object_t> make(arguments_t... args) noexcept
        {
            return make_impl<object_t, false>::make(args...);
        }

        template <typename object_t, typename... arguments_t>
        static inline typename std::enable_if<!std::is_floating_point<object_t>::value, std::shared_ptr<object_t>>::type
            make_shared(arguments_t... args) noexcept
        {
            return make_shared(make_impl<object_t, false>::make(args...));
        }

        template <typename floatp_t, typename kernel_t, typename... arguments_t>
        static inline typename std::enable_if<std::is_floating_point<floatp_t>::value, std::shared_ptr<geometric_object<floatp_t>>>::type
            make_shared(geometry_t geometry, arguments_t... args) noexcept
        {
            // Switch between POINT, LINE, TRIANGLE, TETRAHEDRON, PLANE, RECTANGLE, CUBOID, POLYGON, POLYHEDRON
            switch (static_cast<geometry_t>(geometry))
            {
            case geometry_t::POINT:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<point<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::LINE:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<line<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::TRIANGLE:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<triangle<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::TETRAHEDRON:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<tetrahedron<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::PLANE:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<plane<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::RECTANGLE:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<rectangle<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::CUBOID:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<cuboid<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::POLYGON:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<polygon<floatp_t, kernel_t>, true>::make(args...)));
                break;
            case geometry_t::POLYHEDRON:
                return std::dynamic_pointer_cast<geometric_object<floatp_t>>(make_shared(make_impl<polyhedron<floatp_t, kernel_t>, true>::make(args...)));
                break;
            }

            return nullptr;
        }

        template <typename floatp_t, typename kernel_t>
        static inline std::shared_ptr<geometric_object<floatp_t>> deserialize(const std::vector<char>& serialized)
        {
            const int indicator = static_cast<int>(serialized[0]);

            // Switch between POINT, LINE, TRIANGLE, TETRAHEDRON, PLANE, RECTANGLE, CUBOID, POLYGON, POLYHEDRON
            switch (static_cast<geometry_t>(indicator))
            {
            case geometry_t::POINT:
                return point<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::LINE:
                return line<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::TRIANGLE:
                return triangle<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::TETRAHEDRON:
                return tetrahedron<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::PLANE:
                return plane<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::RECTANGLE:
                return rectangle<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::CUBOID:
                return cuboid<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::POLYGON:
                return polygon<floatp_t, kernel_t>::deserialize(serialized);
                break;
            case geometry_t::POLYHEDRON:
                return polyhedron<floatp_t, kernel_t>::deserialize(serialized);
            }

            throw std::runtime_error(__tpf_error_message("Unknown identifier in serialized geometric object."));
        }
    }
}