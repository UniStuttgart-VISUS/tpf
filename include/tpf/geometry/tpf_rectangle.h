#pragma once

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include <CGAL/Iso_rectangle_2.h>

#include <array>
#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a rectangle in 2d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class rectangle : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = typename kernel_t::Iso_rectangle_2;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point_min">Minimum corner of the rectangle</param>
            /// <param name="point_max">Maximum corner of the rectangle</param>
            /// <param name="render_as_triangles">Render as triangles? (default: polygon)</param>
            rectangle(const point<floatp_t, kernel_t>& point_min, const point<floatp_t, kernel_t>& point_max, bool render_as_triangles = false) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Rectangle to copy from</param>
            rectangle(const rectangle<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Rectangle to copy from</param>
            rectangle(const typename kernel_t::Iso_rectangle_2& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            rectangle<floatp_t, kernel_t>& operator=(const rectangle<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            rectangle<floatp_t, kernel_t>& operator=(const typename kernel_t::Iso_rectangle_2& copy) noexcept;

            /// <summary>
            /// Comparison with another rectangle
            /// </summary>
            /// <param name="other">Other rectangle</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const rectangle<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another rectangle
            /// </summary>
            /// <param name="other">Other rectangle</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const typename kernel_t::Iso_rectangle_2& other) const noexcept;

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
            /// Calculate the area of the rectangle
            /// </summary>
            /// <return>Area of the rectangle</return>
            geometric_float<floatp_t, typename kernel_t::FT> calculate_area() const;

            /// <summary>
            /// Return the minimum point (bottom left)
            /// </summary>
            /// <returns>Minimum point</returns>
            point<floatp_t, kernel_t> get_min_point() const;

            /// <summary>
            /// Return the maximum point (top right)
            /// </summary>
            /// <returns>Maximum point</returns>
            point<floatp_t, kernel_t> get_max_point() const;

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
            /// Return internal representation
            /// </summary>
            /// <returns>Internal representation</returns>
            const typename kernel_t::Iso_rectangle_2& get_internal() const;
            operator const typename kernel_t::Iso_rectangle_2&() const;

        private:
            /// Rectangle
            typename kernel_t::Iso_rectangle_2 _rectangle;

            /// Render as triangle instead as a polygon
            bool render_as_triangles;
        };
    }
}

#include "tpf_rectangle.inl"