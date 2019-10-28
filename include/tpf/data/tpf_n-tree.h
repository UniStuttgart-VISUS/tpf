#pragma once

#include "tpf_position.h"
#include "tpf_tree.h"

#include "../geometry/tpf_point.h"

#include "../policies/tpf_interpolatable.h"

#include <array>
#include <list>
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
        /// <template name="geometry_t">Geometric type used for cells</template>
        /// <template name="num_children">Number of children per node</template>
        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        class n_tree :
            public tree<std::pair<std::shared_ptr<geometry_t>, value_t>>,
            public policies::interpolatable<value_t, geometry::point<float_t, typename geometry_t::kernel_type>>
        {
        public:
            using float_type = float_t;
            using store_type = std::pair<std::shared_ptr<geometry_t>, value_t>;
            using value_type = value_t;
            using kernel_type = typename geometry_t::kernel_type;

            /// <summary>
            /// Node of a tree, storing a value and children depending on derived type
            /// </summary>
            struct node : public tree<store_type>::node
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
                virtual std::shared_ptr<node> set_child(position_t position, const value_t& value) = 0;

                /// <summary>
                /// Set all children of this node
                /// </summary>
                /// <param name="values">New values</param>
                virtual void set_children(const std::array<value_t, num_children>& values) = 0;

                /// <summary>
                /// Get a specific child of this node
                /// </summary>
                /// <param name="position">Position of the child to return</param>
                /// <returns>Requested child</returns>
                std::shared_ptr<node> get_child(position_t position) const;

                /// <summary>
                /// Does this node have the specified child?
                /// </summary>
                /// <returns>True if it contains this child; false otherwise</returns>
                bool has_child(position_t position) const;

                /// <summary>
                /// Get all children attached to this node
                /// </summary>
                /// <returns>All children</returns>
                virtual std::vector<std::shared_ptr<typename tree<store_type>::node>> get_children() const override;

                /// <summary>
                /// Return the number of children
                /// </summary>
                /// <returns>Number of children</returns>
                virtual size_t get_num_children() const override;

            protected:
                // Children, depending on derived type
                std::array<std::shared_ptr<node>, num_children> children;
            };

            /// <summary>
            /// Constructor
            /// </summary>
            n_tree();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="root">Root node</param>
            n_tree(std::shared_ptr<node> root);

            /// <summary>
            /// Return the root of this tree
            /// </summary>
            /// <param name="root">New root node</param>
            /// <returns>Root node</returns>
            std::shared_ptr<node> set_root(std::shared_ptr<node> root);

            /// <summary>
            /// Return the root of this tree
            /// </summary>
            /// <returns>Root node</returns>
            std::shared_ptr<node> get_root() const;

            /// <summary>
            /// Insert nodes
            /// </summary>
            /// <template name="iterator_t">Iterator type for pointing at elements of type pair with [position path, value]</template>
            /// <param name="begin">Iterator pointing at the first element</param>
            /// <param name="end">Iterator pointing past the last element</param>
            template <typename iterator_t>
            void insert_nodes(iterator_t begin, iterator_t end);

            /// <summary>
            /// Find leaf node containing a point
            /// </summary>
            /// <param name="point">Point for which to find a leaf node</param>
            /// <returns>The corresponding leaf node</returns>
            virtual std::pair<std::shared_ptr<node>, std::list<position_t>> find_node(const geometry::point<float_t, kernel_type>& point) const = 0;

            /// <summary>
            /// Interpolate the value at a given point
            /// </summary>
            /// <param name="point">Point at which the value should be interpolated</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate(const geometry::point<float_t, kernel_type>& point) const override;

        protected:
            /// <summary>
            /// Interpolate the values of the given nodes at the center
            /// </summary>
            /// <param name="nodes">Nodes between which interpolation is performed</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate_at_center(const std::vector<std::shared_ptr<typename tree<store_type>::node>>& nodes) const = 0;

            /// <summary>
            /// Interpolate the values at a given point
            /// </summary>
            /// <param name="point">Point at which the value should be interpolated</param>
            /// <param name="nodes">Nodes between which interpolation is performed</param>
            /// <returns>Interpolated value</returns>
            virtual value_t interpolate_at_point(const geometry::point<float_t, kernel_type>& point,
                const std::vector<std::pair<geometry::point<float_t, kernel_type>, value_t>>& nodes) const = 0;

            /// <summary>
            /// Get directions in which interpolation has to be performed
            /// </summary>
            /// <param name="relative_point">Point coordinates relative to the mid point of its cell</param>
            /// <returns>Interpolation directions</returns>
            virtual std::array<neighbor_t, num_children - 1> get_directions(const geometry::point<float_t, kernel_type>& relative_point) const = 0;

        private:
            /// <summary>
            /// Insert a child at a specified position and create all intermediate nodes
            /// </summary>
            /// <param name="position">Position of new child</param>
            /// <param name="value">Value of new child</param>
            /// <returns>The input child</returns>
            std::shared_ptr<node> insert_node(const std::vector<position_t>& position, value_t& value);

            /// <summary>
            /// Get node specified by the position, or nullptr if it does not exist
            /// </summary>
            /// <param name="position">Position of the requested node</param>
            /// <returns>Node corresponding to position, or nullptr if it does not exist</returns>
            std::shared_ptr<node> get_node(const std::vector<position_t>& position) const;
        };
    }
}

#include "tpf_n-tree.inl"