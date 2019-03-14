#pragma once

#include "tpf_tree.h"

#include <memory>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Binary tree data structure
        /// </summary>
        /// <template name="value_t">Storage type of the nodes</template>
        template <typename value_t>
        class binary_tree : public tree<value_t>
        {
        public:
            /// <summary>
            /// Node of a tree, storing a value and its (max. of two) children
            /// </summary>
            struct node : public tree<value_t>::node
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="value">Value to store at this node</param>
                node(const value_t& value);
                node(value_t&& value);

                /// <summary>
                /// Destructor
                /// </summary>
                virtual ~node() override;
                
                /// <summary>
                /// Set the left child of this node
                /// </summary>
                /// <param name="left_child">New left child</param>
                std::shared_ptr<node> set_left_child(std::shared_ptr<node> left_child);

                /// <summary>
                /// Set the right child of this node
                /// </summary>
                /// <param name="right_child">New right child</param>
                std::shared_ptr<node> set_right_child(std::shared_ptr<node> right_child);

                /// <summary>
                /// Set both children of this node
                /// </summary>
                /// <param name="left_child">New left child</param>
                /// <param name="right_child">New right child</param>
                void set_children(std::shared_ptr<node> left_child, std::shared_ptr<node> right_child);

                /// <summary>
                /// Get the left child of this node
                /// </summary>
                /// <returns>Left child</returns>
                std::shared_ptr<node> get_left_child() const;

                /// <summary>
                /// Get the right child of this node
                /// </summary>
                /// <returns>Right child</returns>
                std::shared_ptr<node> get_right_child() const;

                /// <summary>
                /// Does this node have a left child?
                /// </summary>
                /// <returns>True if there is a left child; false otherwise</returns>
                bool has_left_child() const;

                /// <summary>
                /// Does this node right a left child?
                /// </summary>
                /// <returns>True if there is a right child; false otherwise</returns>
                bool has_right_child() const;

                /// <summary>
                /// Get all children attached to this node
                /// </summary>
                /// <returns>All children</returns>
                virtual std::vector<std::shared_ptr<typename tree<value_t>::node>> get_children() const override;

                /// <summary>
                /// Return the number of children
                /// </summary>
                /// <returns>Number of children</returns>
                virtual size_t get_num_children() const override;

            private:
                /// Left and right child
                std::shared_ptr<node> left_child, right_child;
            };

            /// <summary>
            /// Constructor
            /// </summary>
            binary_tree();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="root">Root node</param>
            binary_tree(std::shared_ptr<node> root);

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
        };
    }
}

#include "tpf_binary_tree.inl"