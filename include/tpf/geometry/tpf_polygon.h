#pragma once

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_triangle.h"

#include "../math/tpf_transformer.h"

#include <CGAL/Polygon_2.h>

#include "Eigen/Dense"

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a polygon in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class polygon : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = CGAL::Polygon_2<kernel_t>;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="points">Points of the polygon</param>
            /// <param name="make_convex">Make polygon convex by reordering its vertices</param>
            /// <param name="render_as_triangles">Render as triangles? (default: polygon)</param>
            polygon(const std::vector<point<floatp_t, kernel_t>>& points, bool make_convex, bool render_as_triangles = false) noexcept(false);

            /// <summary>
            /// Constructor with normal
            /// </summary>
            /// <param name="points">Points of the polygon</param>
            /// <param name="normal">Normal of the polygon</param>
            /// <param name="make_convex">Make polygon convex by reordering its vertices</param>
            /// <param name="render_as_triangles">Render as triangles? (default: polygon)</param>
            polygon(std::vector<point<floatp_t, kernel_t>> points, const Eigen::Matrix<floatp_t, 3, 1>& normal, bool make_convex, bool render_as_triangles = false) noexcept(false);

            /// <summary>
            /// Construct a polygon with three vertices
            /// </summary>
            /// <param name="triangle">Triangle</param>
            polygon(const triangle<floatp_t, kernel_t>& triangle) noexcept;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="polygon">Polygon</param>
            polygon(const CGAL::Polygon_2<kernel_t>& polygon) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Polygon to copy from</param>
            polygon(const polygon<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="triangle">Triangle</param>
            /// <returns>This</returns>
            polygon<floatp_t, kernel_t>& operator=(const triangle<floatp_t, kernel_t>& triangle) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="polygon">Polygon</param>
            /// <returns>This</returns>
            polygon<floatp_t, kernel_t>& operator=(const CGAL::Polygon_2<kernel_t>& polygon) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            polygon<floatp_t, kernel_t>& operator=(const polygon<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Comparison with another polygon
            /// </summary>
            /// <param name="other">Other polygon</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const polygon<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another CGAL polygon
            /// </summary>
            /// <param name="other">Other CGAL polygon</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const CGAL::Polygon_2<kernel_t>& other) const noexcept;

            /// <summary>
            /// Clone object (deep copy).
            /// Note that transformations are not supported.
            /// </summary>
            /// <returns>Deep copy</returns>
            virtual std::shared_ptr<geometric_object<floatp_t>> clone(const math::transformer<floatp_t, 3>& trafo = math::transformer<floatp_t, 3>::unit()) const;

            /// <summary>
            /// Note that transformations are not supported
            /// </summary>
            /// <param name="trafo">Transformer</param>
            virtual geometric_object<floatp_t>& transform(const math::transformer<floatp_t, 3>& trafo);

            /// <summary>
            /// Calculate the polygon area
            /// </summary>
            /// <returns>Polygon area</returns>
            geometric_float<floatp_t, typename kernel_t::FT> calculate_area() const;

            /// <summary>
            /// Calculate the centroid of the polygon, if it is convex
            /// </summary>
            /// <returns>Centroid</returns>
            point<floatp_t, kernel_t> calculate_centroid() const;

            /// <summary>
            /// Answer if the polygon is convex or not
            /// </summary>
            /// <returns>True if convex; false otherwise</returns>
            bool is_convex() const;

            /// <summary>
            /// Return number of points
            /// </summary>
            /// <returns>Number of points</returns>
            virtual std::size_t get_num_points() const;

            /// <summary>
            /// Return points
            /// </summary>
            /// <returns>Points</returns>
            virtual std::vector<Eigen::Matrix<floatp_t, 3, 1>> get_points() const;

            /// <summary>
            /// Return number of cells
            /// </summary>
            /// <returns>Number of cells</returns>
            virtual std::size_t get_num_cells() const;

            /// <summary>
            /// Return cells
            /// </summary>
            /// <returns>Cells</returns>
            virtual std::vector<std::vector<std::size_t>> get_cells() const;

            /// <summary>
            /// Return the centroid or center of mass of the object
            /// </summary>
            /// <returns>Centroid</returns>
            virtual Eigen::Matrix<floatp_t, 3, 1> get_centroid() const;

            /// <summary>
            /// Serialize object to a binary representation
            /// </summary>
            /// <returns>Binary representation</returns>
            virtual std::vector<char> serialize() const;

            /// <summary>
            /// Deserialize object from a binary representation
            /// </summary>
            /// <param name="serialized">Serialized object</param>
            /// <returns>Deserialized object</returns>
            static std::shared_ptr<geometric_object<floatp_t>> deserialize(const std::vector<char>& serialized);

            /// <summary>
            /// Answer the type of this geometric object
            /// </summary>
            /// <returns>Type of this geometric object</returns>
            virtual geometry_t get_type() const;

            /// <summary>
            /// Return internal representation
            /// </summary>
            /// <returns>Internal representation</returns>
            const CGAL::Polygon_2<kernel_t>& get_internal() const;
            operator const CGAL::Polygon_2<kernel_t>&() const;

            /// <summary>
            /// Get transformation
            /// </summary>
            /// <returns>Transformation</returns>
            const math::transformer<floatp_t, 3>& get_transformation() const;

        private:
            /// <summary>
            /// Calculate normal from points
            /// </summary>
            /// <param name="points">Input points</param>
            /// <returns>Normal</returns>
            Eigen::Matrix<floatp_t, 3, 1> calculate_normal(const std::vector<point<floatp_t, kernel_t>>& points) const noexcept(false);

            /// Polygon
            CGAL::Polygon_2<kernel_t> _polygon;

            /// Transformation for support of polygons in 3D space
            math::transformer<floatp_t, 3> transformation;

            /// Render as triangles?
            bool render_as_triangles;
        };
    }
}

#include "tpf_polygon.inl"
