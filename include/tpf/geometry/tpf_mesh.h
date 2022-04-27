#pragma once

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_triangle.h"

#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <initializer_list>
#include <memory>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Represents a mesh in 3d space
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        class mesh : public geometric_object<floatp_t>
        {
        public:
            using kernel_type = kernel_t;
            using internal_type = CGAL::Surface_mesh<typename kernel_t::Point_3>;

            /// <summary>
            /// Create empty mesh
            /// </summary>
            explicit mesh(bool repair = false) noexcept;

            /// <summary>
            /// Construct from existing mesh
            /// </summary>
            /// <param name="mesh">Mesh</param>
            mesh(const CGAL::Surface_mesh<typename kernel_t::Point_3>& mesh, bool repair = false) noexcept;

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Mesh to copy from</param>
            mesh(const mesh<floatp_t, kernel_t>& copy) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="mesh">Mesh</param>
            /// <returns>This</returns>
            mesh<floatp_t, kernel_t>& operator=(const CGAL::Surface_mesh<typename kernel_t::Point_3>& mesh) noexcept;

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            mesh<floatp_t, kernel_t>& operator=(const mesh<floatp_t, kernel_t>& copy) noexcept;

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
            /// Add a point, which is not connected until it is used in a face
            /// </summary>
            /// <param name="point">New point to add</param>
            /// <returns>Index of the point</returns>
            CGAL::SM_Vertex_index add_point(const point<floatp_t, kernel_t>& point);

            /// <summary>
            /// Add a face by connecting three points
            /// (throws an exception if the point index is invalid)
            /// </summary>
            /// <param name="v0">Point index or new point to add</param>
            /// <param name="v1">Point index or new point to add</param>
            /// <param name="v2">Point index or new point to add</param>
            void add_face(CGAL::SM_Vertex_index v0, CGAL::SM_Vertex_index v1, CGAL::SM_Vertex_index v2);
            void add_face(const point<floatp_t, kernel_t>& v0, const point<floatp_t, kernel_t>& v1, const point<floatp_t, kernel_t>& v2);

            /// <summary>
            /// Add a face by connecting multiple co-planar points
            /// (throws an exception if points are not co-planar or number of points is less than three)
            /// </summary>
            /// <param name="v">Point indices or new points to add</param>
            void add_face(const std::vector<CGAL::SM_Vertex_index>& v);
            void add_face(const std::initializer_list<CGAL::SM_Vertex_index>& v);
            void add_face(const std::vector<point<floatp_t, kernel_t>>& v);
            void add_face(const std::initializer_list<point<floatp_t, kernel_t>>& v);

            /// <summary>
            /// Add a triangle face
            /// </summary>
            /// <param name="triangle">Triangle to add</param>
            void add_face(const triangle<floatp_t, kernel_t>& triangle);

            /// <summary>
            /// Merge this mesh with another, resulting in a single mesh
            /// </summary>
            /// <param name="other"></param>
            std::shared_ptr<mesh<floatp_t, kernel_t>> merge(const mesh<floatp_t, kernel_t>& other) const;

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
            /// Return the centroid or center of mass of the object
            /// </summary>
            /// <returns>Centroid</returns>
            virtual Eigen::Matrix<floatp_t, 3, 1> get_centroid() const;

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
            const CGAL::Surface_mesh<typename kernel_t::Point_3>& get_internal() const;
            operator const CGAL::Surface_mesh<typename kernel_t::Point_3>&() const;

            /// <summary>
            /// Repair the mesh if needed
            /// </summary>
            void repair() const;

        private:
            /// Triangle
            CGAL::Surface_mesh<typename kernel_t::Point_3> _mesh;

            /// Should the mesh be repaired if necessary?
            bool perform_repair;

            /// Flag signalling need for repair
            bool valid;
        };
    }
}

#include "tpf_mesh.inl"
