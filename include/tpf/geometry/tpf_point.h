#pragma once

#include "tpf_geometric_object.h"

#include <CGAL/Point_3.h>

#include "Eigen/Dense"

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a point in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class point : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = typename kernel_t::Point_3;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="x">X value</param>
            /// <param name="y">Y value</param>
            /// <param name="z">Z value</param>
            point(floatp_t x, floatp_t y, floatp_t z) noexcept;
            
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point">Point</param>
            point(const Eigen::Matrix<floatp_t, 3, 1>& point) noexcept;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point">Point</param>
            point(const typename kernel_t::Point_3& point) noexcept;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point">Point</param>
            point(const typename kernel_t::Point_2& point, floatp_t z = static_cast<floatp_t>(0.0)) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Point to copy from</param>
            point(const point<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="point">Point</param>
            /// <returns>This</returns>
            point<floatp_t, kernel_t>& operator=(const typename kernel_t::Point_3& point) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="point">Point</param>
            /// <returns>This</returns>
            point<floatp_t, kernel_t>& operator=(const typename kernel_t::Point_2& point) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            point<floatp_t, kernel_t>& operator=(const point<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Assignment from a vector
            /// </summary>
            /// <param name="point">Point</param>
            /// <returns>This</returns>
            point<floatp_t, kernel_t>& operator=(const Eigen::Matrix<floatp_t, 3, 1>& point) noexcept;

            /// <summary>
            /// Comparison with another point
            /// </summary>
            /// <param name="other">Other point</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const point<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another CGAL triangle
            /// </summary>
            /// <param name="other">Other CGAL triangle</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const typename kernel_t::Point_3& other) const noexcept;

            /// <summary>
            /// Clone object (deep copy)
            /// </summary>
            /// <returns>Deep copy</returns>
            virtual std::shared_ptr<geometric_object<floatp_t>> clone(const math::transformer<floatp_t, 3>& trafo = math::transformer<floatp_t, 3>::unit()) const override;

            /// <summary>
            /// Transform this object using a transformer
            /// </summary>
            /// <param name="trafo">Transformer</param>
            virtual geometric_object<floatp_t>& transform(const math::transformer<floatp_t, 3>& trafo) override;

            /// <summary>
            /// Get vertex
            /// </summary>
            /// <returns>Vertex</returns>
            Eigen::Matrix<floatp_t, 3, 1> get_vertex() const;
            operator Eigen::Matrix<floatp_t, 3, 1>() const;

            /// <summary>
            /// Get vector element at the given position
            /// </summary>
            /// <param name="i">Element position</param>
            /// <returns>Vector element</returns>
            floatp_t operator[](std::size_t i) const;

            /// <summary>
            /// Return number of points
            /// </summary>
            /// <returns>Number of points</returns>
            virtual std::size_t get_num_points() const override;

            /// <summary>
            /// Return points
            /// </summary>
            /// <returns>Points</returns>
            virtual std::vector<Eigen::Matrix<floatp_t, 3, 1>> get_points() const override;

            /// <summary>
            /// Return number of cells
            /// </summary>
            /// <returns>Number of cells</returns>
            virtual std::size_t get_num_cells() const override;

            /// <summary>
            /// Return cells
            /// </summary>
            /// <returns>Cells</returns>
            virtual std::vector<std::vector<std::size_t>> get_cells() const override;

            /// <summary>
            /// Serialize object to a binary representation
            /// </summary>
            /// <returns>Binary representation</returns>
            virtual std::vector<char> serialize() const override;

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
            const typename kernel_t::Point_3& get_internal() const;
            operator const typename kernel_t::Point_3&() const;

            explicit operator typename kernel_t::Point_2() const;

        private:
            /// Point
            typename kernel_t::Point_3 _point;
        };
    }
}

#include "tpf_point.inl"