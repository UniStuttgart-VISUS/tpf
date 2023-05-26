#pragma once

#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a line strip in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class line_strip : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = std::vector<point<floatp_t, kernel_t>>;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point_1">First point on the line</param>
            /// <param name="point_2">Second point on the line</param>
            line_strip(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2) noexcept;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="points">Points in correct order, describing a line strip</param>
            line_strip(const std::vector<point<floatp_t, kernel_t>>& points) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Line to copy from</param>
            line_strip(const line_strip<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            line_strip<floatp_t, kernel_t>& operator=(const line_strip<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Comparison with another line strip
            /// </summary>
            /// <param name="other">Other line strip</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const line_strip<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Clone object (deep copy)
            /// </summary>
            /// <returns>Deep copy</returns>
            virtual std::shared_ptr<geometric_object<floatp_t>> clone(const math::transformer<floatp_t, 3>& trafo = math::transformer<floatp_t, 3>::unit()) const;

            /// <summary>
            /// Transform this object using a transformer
            /// </summary>
            /// <param name="trafo">Transformer</param>
            virtual geometric_object<floatp_t>& transform(const math::transformer<floatp_t, 3>& trafo);

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
            const std::vector<point<floatp_t, kernel_t>>& get_internal() const;
            operator const std::vector<point<floatp_t, kernel_t>>&() const;

        private:
            /// Line strip
            std::vector<point<floatp_t, kernel_t>> _line_strip;
        };
    }
}

#include "tpf_line_strip.inl"
