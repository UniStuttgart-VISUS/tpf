#include "tpf_mesh.h"

#include "tpf_geometric_float.h"
#include "tpf_geometric_object.h"
#include "tpf_point.h"
#include "tpf_triangle.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../log/tpf_log.h"

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>::mesh() noexcept : valid(true)
        {}

        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>::mesh(const CGAL::Surface_mesh<typename kernel_t::Point_3>& mesh) noexcept : valid(false)
        {
            this->_mesh = mesh;
        }

        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>::mesh(const mesh<floatp_t, kernel_t>& copy) noexcept
        {
            this->_mesh = copy._mesh;
            this->valid = copy.valid;
        }

        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>& mesh<floatp_t, kernel_t>::operator=(const CGAL::Surface_mesh<typename kernel_t::Point_3>& mesh) noexcept
        {
            geometric_object<floatp_t>::operator=(mesh);

            this->_mesh = mesh;
            this->valid = false;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>& mesh<floatp_t, kernel_t>::operator=(const mesh<floatp_t, kernel_t>& copy) noexcept
        {
            geometric_object<floatp_t>::operator=(copy);

            this->_mesh = copy._mesh;
            this->valid = copy.valid;

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<geometric_object<floatp_t>> mesh<floatp_t, kernel_t>::clone(const math::transformer<floatp_t, 3>& trafo) const
        {
            auto copy = std::make_shared<mesh<floatp_t, kernel_t>>(*this);
            copy->transform(trafo);

            return copy;
        }

        template <typename floatp_t, typename kernel_t>
        inline geometric_object<floatp_t>& mesh<floatp_t, kernel_t>::transform(const math::transformer<floatp_t, 3>& trafo)
        {
            if (!trafo.is_unit())
            {
                const auto vertices = this->_mesh.vertices();

                for (auto it = vertices.begin(); it != vertices.end(); ++it)
                {
                    auto& p = this->_mesh.point(*it);

                    const auto original = point<floatp_t, kernel_t>(p).get_vertex();
                    const auto transformation = trafo.transform(original) - original;

                    typename kernel_t::Vector_3 vector(transformation[0], transformation[1], transformation[2]);

                    p += vector;
                }
            }

            return *this;
        }

        template <typename floatp_t, typename kernel_t>
        inline CGAL::SM_Vertex_index mesh<floatp_t, kernel_t>::add_point(const point<floatp_t, kernel_t>& point)
        {
            this->valid = false;

            return this->_mesh.add_vertex(point.get_internal());
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(CGAL::SM_Vertex_index v0, CGAL::SM_Vertex_index v1, CGAL::SM_Vertex_index v2)
        {
            const auto face_index = this->_mesh.add_face(v0, v1, v2);

            if (face_index == typename CGAL::Surface_mesh<typename kernel_t::Point_3>::null_face())
            {
                throw std::runtime_error(__tpf_error_message("Face could not be added to the mesh."));
            }

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const point<floatp_t, kernel_t>& p0, const point<floatp_t, kernel_t>& p1, const point<floatp_t, kernel_t>& p2)
        {
            const auto v0 = add_point(p0);
            const auto v1 = add_point(p1);
            const auto v2 = add_point(p2);

            add_face(v0, v1, v2);

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const std::vector<CGAL::SM_Vertex_index>& indices)
        {
            const auto face_index = this->_mesh.add_face(indices);

            if (face_index == typename CGAL::Surface_mesh<typename kernel_t::Point_3>::null_face())
            {
                throw std::runtime_error(__tpf_error_message("Face could not be added to the mesh."));
            }

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const std::initializer_list<CGAL::SM_Vertex_index>& indices)
        {
            const auto face_index = this->_mesh.add_face(indices);

            if (face_index == typename CGAL::Surface_mesh<typename kernel_t::Point_3>::null_face())
            {
                throw std::runtime_error(__tpf_error_message("Face could not be added to the mesh."));
            }

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const std::vector<point<floatp_t, kernel_t>>& points)
        {
            std::vector<CGAL::SM_Vertex_index> indices;
            indices.reserve(v.size());

            for (const auto& p : points)
            {
                indices.push_back(add_point(p));
            }

            add_face(indices);

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const std::initializer_list<point<floatp_t, kernel_t>>& points)
        {
            std::vector<CGAL::SM_Vertex_index> indices;
            indices.reserve(v.size());

            for (const auto& p : points)
            {
                indices.push_back(add_point(p));
            }

            add_face(indices);

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::add_face(const triangle<floatp_t, kernel_t>& triangle)
        {
            add_face(
                point<floatp_t, kernel_t>(triangle.get_internal().vertex(0)),
                point<floatp_t, kernel_t>(triangle.get_internal().vertex(1)),
                point<floatp_t, kernel_t>(triangle.get_internal().vertex(2)));

            this->valid = false;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::shared_ptr<mesh<floatp_t, kernel_t>> mesh<floatp_t, kernel_t>::merge(const mesh<floatp_t, kernel_t>& other) const
        {
            auto merged = std::static_pointer_cast<mesh<floatp_t, kernel_t>>(other.clone());

            std::map<CGAL::SM_Vertex_index, CGAL::SM_Vertex_index> index_map;

            // Insert own points and store their indices in a map
            const auto vertices = this->_mesh.vertices();

            for (auto it = vertices.begin(); it != vertices.end(); ++it)
            {
                const auto index = *it;

                index_map[index] = merged->add_point(point<floatp_t, kernel_t>(this->_mesh.point(*it)));
            }

            // Add own faces using the index map
            const auto faces = this->_mesh.faces();

            for (auto face_it = faces.begin(); face_it != faces.end(); ++face_it)
            {
                const auto vertices = this->_mesh.vertices_around_face(this->_mesh.halfedge(*face_it));
                std::vector<CGAL::SM_Vertex_index> new_ids;

                for (auto vertex_it = vertices.begin(); vertex_it != vertices.end(); ++vertex_it)
                {
                    new_ids.push_back(index_map.at(*vertex_it));
                }

                merged->add_face(new_ids);
            }

            return merged;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t mesh<floatp_t, kernel_t>::get_num_points() const
        {
            repair();

            return this->_mesh.number_of_vertices();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<Eigen::Matrix<floatp_t, 3, 1>> mesh<floatp_t, kernel_t>::get_points() const
        {
            repair();

            std::vector<Eigen::Matrix<floatp_t, 3, 1>> points;
            points.reserve(get_num_points());

            const auto vertices = this->_mesh.vertices();

            for (auto it = vertices.begin(); it != vertices.end(); ++it)
            {
                points.push_back(point<floatp_t, kernel_t>(this->_mesh.point(*it)).get_vertex());
            }

            return points;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::size_t mesh<floatp_t, kernel_t>::get_num_cells() const
        {
            repair();

            return this->_mesh.number_of_faces();
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<std::vector<std::size_t>> mesh<floatp_t, kernel_t>::get_cells() const
        {
            repair();

            std::vector<std::vector<std::size_t>> cells;
            cells.reserve(get_num_cells());

            const auto faces = this->_mesh.faces();

            for (auto face_it = faces.begin(); face_it != faces.end(); ++face_it)
            {
                const auto vertices = this->_mesh.vertices_around_face(this->_mesh.halfedge(*face_it));

                cells.emplace_back();

                for (auto vertex_it = vertices.begin(); vertex_it != vertices.end(); ++vertex_it)
                {
                    cells.back().push_back(*vertex_it);
                }
            }

            return cells;
        }

        template <typename floatp_t, typename kernel_t>
        inline std::vector<char> mesh<floatp_t, kernel_t>::serialize() const
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        std::shared_ptr<geometric_object<floatp_t>> mesh<floatp_t, kernel_t>::deserialize(const std::vector<char>& serialized)
        {
            throw exception::not_implemented_exception();
        }

        template <typename floatp_t, typename kernel_t>
        inline const CGAL::Surface_mesh<typename kernel_t::Point_3>& mesh<floatp_t, kernel_t>::get_internal() const
        {
            return this->_mesh;
        }

        template <typename floatp_t, typename kernel_t>
        inline mesh<floatp_t, kernel_t>::operator const CGAL::Surface_mesh<typename kernel_t::Point_3>& () const
        {
            return this->_mesh;
        }

        template <typename floatp_t, typename kernel_t>
        inline void mesh<floatp_t, kernel_t>::repair() const
        {
            if (!this->valid)
            {
                CGAL::Polygon_mesh_processing::remove_isolated_vertices(const_cast<CGAL::Surface_mesh<typename kernel_t::Point_3>&>(this->_mesh));
                CGAL::Polygon_mesh_processing::stitch_borders(const_cast<CGAL::Surface_mesh<typename kernel_t::Point_3>&>(this->_mesh));
            }

            const_cast<bool&>(this->valid) = true;
        }
    }
}
