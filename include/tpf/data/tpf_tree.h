#pragma once

#include <list>
#include <memory>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Tree data structure
        /// </summary>
        /// <template name="value_t">Storage type of the nodes</template>
        template <typename value_t>
        class tree
        {
        public:
            using value_type = value_t;

            /// <summary>
            /// Abstract node of a tree
            /// </summary>
            struct node
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                node();

                /// <summary>
                /// Destructor
                /// </summary>
                virtual ~node() = 0;

                /// <summary>
                /// Is this node a leaf?
                /// </summary>
                /// <returns>True if it is a leaf (has no children); false otherwise</returns>
                bool is_leaf() const;

                /// <summary>
                /// Access the value stored at this node
                /// </summary>
                /// <returns></returns>
                const value_t& get_value() const;
                value_t& get_value();

                /// <summary>
                /// Get all children attached to this node
                /// </summary>
                /// <returns>All children</returns>
                virtual std::vector<std::shared_ptr<node>> get_children() const = 0;

                /// <summary>
                /// Return the number of children
                /// </summary>
                /// <returns>Number of children</returns>
                virtual size_t get_num_children() const = 0;

                /// <summary>
                /// Get the parent of this node
                /// </summary>
                /// <returns>Parent node</returns>
                node* get_parent() const;

            protected:
                /// <summary>
                /// Set the parent of this node
                /// </summary>
                /// <param name="parent">New parent node</param>
                void set_parent(node* parent);

            private:
                /// Value stored at the node
                value_t value;

                /// Parent node
                node* parent;
            };

            /// <summary>
            /// Constructor
            /// </summary>
            tree();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="root">Root node</param>
            tree(std::shared_ptr<node> root);

            /// <summary>
            /// Set and return the root of this tree
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
            /// Get all paths from the root to the leaves
            /// </summary>
            /// <returns>Paths from the root to the leaves</returns>
            std::vector<std::list<std::shared_ptr<node>>> get_paths() const;

            /// <summary>
            /// Traverse tree and return all leaf nodes in depth-first-order
            /// </summary>
            /// <returns>List of all leaf nodes</returns>
            std::list<std::shared_ptr<node>> get_leaf_nodes() const;

            /// <summary>
            /// Perform breadth-first-search and return the corresponding list of nodes
            /// </summary>
            /// <returns>List of nodes in breadth-first-search order</returns>
            std::list<std::shared_ptr<node>> breadth_first_search() const;

            /// <summary>
            /// Perform depth-first-search and return the corresponding list of nodes
            /// </summary>
            /// <returns>List of nodes in depth-first-search order</returns>
            std::list<std::shared_ptr<node>> depth_first_search() const;

        private:
            /// <summary>
            /// Get all paths from the given root to the leaves
            /// </summary>
            /// <param name="root">Root of the subtree</param>
            /// <returns>Paths for this subtree</returns>
            std::vector<std::list<std::shared_ptr<node>>> get_paths(std::shared_ptr<node> root) const;

            /// <summary>
            /// Perform depth-first-search and return the corresponding list of nodes
            /// </summary>
            /// <returns>List of nodes in depth-first-search order</returns>
            std::list<std::shared_ptr<node>> depth_first_search(std::shared_ptr<node> root) const;

            /// Root of this tree
            std::shared_ptr<node> root;
        };
    }
}

#include "tpf_tree.inl"