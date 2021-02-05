#pragma once

#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include "Eigen/Dense"

#include <CGAL/Plane_3.h>

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a plane in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class plane : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = typename kernel_t::Plane_3;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point">Up point on the plane</param>
            /// <param name="normal">Plane normal</param>
            plane(const point<floatp_t, kernel_t>& point, const Eigen::Matrix<floatp_t, 3, 1>& normal);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point_1">First point on the plane</param>
            /// <param name="point_2">Second point on the plane</param>
            /// <param name="point_3">Third point on the plane</param>
            plane(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2, const point<floatp_t, kernel_t>& point_3) noexcept;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="plane">Plane</param>
            plane(const typename kernel_t::Plane_3& plane) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Plane to copy from</param>
            plane(const plane<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="plane">Plane</param>
            /// <returns>This</returns>
            plane<floatp_t, kernel_t>& operator=(const typename kernel_t::Plane_3& plane) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            plane<floatp_t, kernel_t>& operator=(const plane<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Comparison with another plane
            /// </summary>
            /// <param name="other">Other plane</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const plane<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another CGAL plane
            /// </summary>
            /// <param name="other">Other CGAL plane</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const typename kernel_t::Plane_3& other) const noexcept;

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
            const typename kernel_t::Plane_3& get_internal() const;
            operator const typename kernel_t::Plane_3&() const;

        private:
            /// Plane
            typename kernel_t::Plane_3 _plane;
        };
    }
}

#include "tpf_plane.inl"