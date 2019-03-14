#pragma once

#include "tpf_binary_tree.h"

#include "../geometry/tpf_geometric_object.h"
#include "../geometry/tpf_plane.h"
#include "../geometry/tpf_point.h"
#include "../geometry/tpf_polygon.h"

#include "Eigen/Dense"

#include <memory>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Binary space partitioning (BSP) tree data structure storing at each node a plane and a normal
        /// </summary>
        /// <template name="float_t">Floating point type of the plane</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename float_t, typename kernel_t = geometry::default_kernel_t>
        class bsp_tree : public binary_tree<std::pair<std::shared_ptr<geometry::plane<float_t, kernel_t>>, Eigen::Matrix<float_t, 3, 1>>>
        {
        public:
            using float_type = float_t;
            using value_type = std::pair<std::shared_ptr<geometry::plane<float_t, kernel_t>>, Eigen::Matrix<float_t, 3, 1>>;
            using input_value_type = std::pair<std::shared_ptr<geometry::polygon<float_t, kernel_t>>, Eigen::Matrix<float_t, 3, 1>>;
            using node_type = typename binary_tree<value_type>::node;
            using kernel_type = kernel_t;

            /// <summary>
            /// Constructor
            /// </summary>
            bsp_tree();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="root">Root node</param>
            bsp_tree(std::shared_ptr<node_type> root);

            /// <summary>
            /// Create the BSP tree from a set of polygons and their normals, denoting the faces of a polyhedron
            /// </summary>
            /// <template name="forward_iterator_t">Forward iterator type</template>
            /// <param name="begin">Iterator pointing to the first element of type pair[shared_ptr[plane], normal]</param>
            /// <param name="end">Iterator pointing past the last element of type pair[shared_ptr[plane], normal]</param>
            template <typename forward_iterator_t>
            void create_bsp_tree(forward_iterator_t begin, forward_iterator_t end);

            /// <summary>
            /// Is this point inside?
            /// </summary>
            /// <param name="point">Point to check</param>
            /// <returns>True if inside; false otherwise</returns>
            bool is_inside(const geometry::point<float_t, kernel_t>& point) const;

            /// <summary>
            /// Extract convex sub trees from this BSP tree
            /// </summary>
            /// <returns>Convex sub trees</returns>
            std::vector<bsp_tree<float_t, kernel_t>> extract_convex_subtrees() const;

        private:
            /// <summary>
            /// Construct an exact partitioning plane from the given points
            /// </summary>
            /// <param name="points">Three input points</param>
            /// <param name="normal">Approximate normal</param>
            /// <returns>Partitioning plane</returns>
            geometry::plane<float_t, kernel_t> construct_partitioning_plane(const std::vector<Eigen::Matrix<float_t, 3, 1>>& points, const Eigen::Matrix<float_t, 3, 1>& normal) const;
        };
    }
}

#include "tpf_bsp_tree.inl"