#pragma once

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_tetrahedron.h"

#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a convex polyhedron in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class polyhedron : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = std::vector<tetrahedron<floatp_t, kernel_t>>;

            /// <summary>
            /// Constructor for a convex polyhedron
            /// </summary>
            /// <param name="points">Edges of the polyhedron</param>
            polyhedron(const std::vector<point<floatp_t, kernel_t>>& points) noexcept(false);

            /// <summary>
            /// Constructor for a convex polyhedron
            /// </summary>
            /// <param name="tetrahedra">Tetrahedra</param>
            polyhedron(const std::vector<tetrahedron<floatp_t, kernel_t>>& tetrahedra) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Polyhedron to copy from</param>
            polyhedron(const polyhedron<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="tetrahedra">Tetrahedra</param>
            /// <returns>This</returns>
            polyhedron<floatp_t, kernel_t>& operator=(const std::vector<tetrahedron<floatp_t, kernel_t>>& tetrahedra) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            polyhedron<floatp_t, kernel_t>& operator=(const polyhedron<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Comparison with another polyhedron
            /// </summary>
            /// <param name="other">Other polyhedron</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const polyhedron<floatp_t, kernel_t>& other) const noexcept;

            /// <summary>
            /// Comparison with another set of tetrahedra
            /// </summary>
            /// <param name="other">Other set of tetrahedra</param>
            /// <returns>True if same; false otherwise</returns>
            bool operator==(const std::vector<tetrahedron<floatp_t, kernel_t>>& other) const noexcept;

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
            /// Calculate the volume of this polyhedron
            /// </summary>
            /// <returns>Volume of this polyhedron</returns>
            geometric_float<floatp_t, typename kernel_t::FT> calculate_volume() const;

            /// <summary>
            /// Calculate the centroid of this polyhedron
            /// </summary>
            /// <returns>Centroid of this polyhedron</returns>
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
            const std::vector<tetrahedron<floatp_t, kernel_t>>& get_internal() const;
            operator const std::vector<tetrahedron<floatp_t, kernel_t>>&() const;

        private:
            /// Tetrahedra
            std::vector<tetrahedron<floatp_t, kernel_t>> _tetrahedra;
        };
    }
}

#include "tpf_polyhedron.inl"