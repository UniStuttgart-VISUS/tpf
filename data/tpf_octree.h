#pragma once

#include "tpf_n-tree.h"
#include "tpf_position.h"

#include "../geometry/tpf_cuboid.h"
#include "../geometry/tpf_geometric_object.h"
#include "../geometry/tpf_point.h"

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Octree data structure
        /// </summary>
        /// <template name="float_t">Floating point type for positions</template>
        /// <template name="value_t">Value type to store at each cell center</template>
        /// <template name="kernel_t">Internal CGAL kernel</template>
        template <typename float_t, typename value_t, typename kernel_t = geometry::default_kernel_t>
        class octree : public n_tree<float_t, value_t, geometry::cuboid<float_t, kernel_t>, 8>
        {
        public:
            using base_type = n_tree<float_t, value_t, geometry::cuboid<float_t, kernel_t>, 8>;
            using float_type = float_t;
            using store_type = std::pair<std::shared_ptr<geometry::cuboid<float_t, kernel_t>>, value_t>;
            using value_type = value_t;
            using kernel_type = kernel_t;

            /// <summary>
            /// Node of a tree, storing a value and up to eight children
            /// </summary>
            struct node : public base_type::node
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="value">Value to store at this node</param>
                node(const store_type& value);
                node(store_type&& value);

                /// <summary>
                /// Set a specific child of this node
                /// </summary>
                /// <param name="position">Position of new child</param>
                /// <param name="value">Value of new child</param>
                virtual std::shared_ptr<typename base_type::node> set_child(position_t position, const value_t& value) override;

                /// <summary>
                /// Set all children of this node
                /// </summary>
                /// <param name="values">New values</param>
                virtual void set_children(const std::array<value_t, 8>& values) override;
            };

            /// <summary>
            /// Constructor
            /// </summary>
            octree();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="root">Root node</param>
            octree(std::shared_ptr<node> root);

            /// <summary>
            /// Find leaf node containing a point
            /// </summary>
            /// <param name="point">Point for which to find a leaf node</param>
            /// <returns>The corresponding leaf node</returns>
            virtual std::pair<std::shared_ptr<typename base_type::node>, std::list<position_t>> find_node(const geometry::point<float_t, kernel_type>& point) const override;

        protected:
            /// <summary>
            /// Interpolate the values of the given nodes at the center
            /// </summary>
            /// <param name="nodes">Nodes between which interpolation is performed</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate_at_center(const std::vector<std::shared_ptr<typename tree<store_type>::node>>& nodes) const override;

            /// <summary>
            /// Interpolate the values at a given point
            /// </summary>
            /// <param name="point">Point at which the value should be interpolated</param>
            /// <param name="nodes">Nodes between which interpolation is performed</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate_at_point(const geometry::point<float_t, kernel_type>& point,
                const std::vector<std::pair<geometry::point<float_t, kernel_t>, value_t>>& nodes) const override;

            /// <summary>
            /// Get directions in which interpolation has to be performed
            /// </summary>
            /// <param name="relative_point">Point coordinates relative to the mid point of its cell</param>
            /// <returns>Interpolation directions</returns>
            virtual std::array<neighbor_t, 7> get_directions(const geometry::point<float_t, kernel_type>& relative_point) const override;
        };
    }
}

#include "tpf_octree.inl"