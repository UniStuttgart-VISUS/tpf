#pragma once

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"

#include <CGAL/Tetrahedron_3.h>

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a tetrahedron in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class tetrahedron : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = typename kernel_t::Tetrahedron_3;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="points">Edges of the tetrahedron</param>
            tetrahedron(const std::vector<point<floatp_t, kernel_t>>& points) noexcept(false);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point_1">Edge of the tetrahedron</param>
            /// <param name="point_2">Edge of the tetrahedron</param>
            /// <param name="point_3">Edge of the tetrahedron</param>
            /// <param name="point_4">Edge of the tetrahedron</param>
            tetrahedron(const point<floatp_t, kernel_t>& point_1, const point<floatp_t, kernel_t>& point_2,
                const point<floatp_t, kernel_t>& point_3, const point<floatp_t, kernel_t>& point_4) noexcept(false);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="tetrahedron">Tetrahedron</param>
            tetrahedron(const typename kernel_t::Tetrahedron_3& tetrahedron) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Tetrahedron to copy from</param>
            tetrahedron(const tetrahedron<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="tetrahedron">Tetrahedron</param>
            /// <returns>This</returns>
            tetrahedron<floatp_t, kernel_t>& operator=(const typename kernel_t::Tetrahedron_3& tetrahedron) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            tetrahedron<floatp_t, kernel_t>& operator=(const tetrahedron<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Comparison with another tetrahedron
            /// </summary>
            /// <param name="other">Other tetrahedron</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const tetrahedron<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another CGAL tetrahedron
            /// </summary>
            /// <param name="other">Other CGAL tetrahedron</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const typename kernel_t::Tetrahedron_3& other) const noexcept;

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
            /// Calculate the volume of this tetrahedron
            /// </summary>
            /// <returns>Volume of this tetrahedron</returns>
            geometric_float<floatp_t, typename kernel_t::FT> calculate_volume() const;

            /// <summary>
            /// Calculate the centroid of this tetrahedron
            /// </summary>
            /// <returns>Centroid</returns>
            point<floatp_t, kernel_t> calculate_centroid() const;

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
            const typename kernel_t::Tetrahedron_3& get_internal() const;
            operator const typename kernel_t::Tetrahedron_3&() const;

        private:
            /// Tetrahedron
            typename kernel_t::Tetrahedron_3 _tetrahedron;
        };
    }
}

#include "tpf_tetrahedron.inl"
