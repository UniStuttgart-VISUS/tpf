#include "tpf_binary_tree.h"

#include "tpf_tree.h"

#include <memory>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename value_t>
        inline binary_tree<value_t>::node::node(const value_t& value)
        {
            tree<value_t>::node::get_value() = value;

            this->left_child = nullptr;
            this->right_child = nullptr;
        }

        template <typename value_t>
        inline binary_tree<value_t>::node::~node()
        {
            if (this->left_child != nullptr)
            {
                this->left_child->set_parent(nullptr);
            }

            if (this->right_child != nullptr)
            {
                this->right_child->set_parent(nullptr);
            }
        }

        template <typename value_t>
        inline binary_tree<value_t>::node::node(value_t&& value)
        {
            tree<value_t>::node::get_value() = std::move(value);

            this->left_child = nullptr;
            this->right_child = nullptr;
        }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::node::set_left_child(std::shared_ptr<node> left_child)
        {
            left_child->set_parent(this);
            return this->left_child = left_child;
        }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::node::set_right_child(std::shared_ptr<node> right_child)
        {
            right_child->set_parent(this);
            return this->right_child = right_child;
        }

        template <typename value_t>
        inline void binary_tree<value_t>::node::set_children(std::shared_ptr<node> left_child, std::shared_ptr<node> right_child)
        {
            left_child->set_parent(this);
            right_child->set_parent(this);

            this->left_child = left_child;
            this->right_child = right_child;
        }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::node::get_left_child() const
        {
            return this->left_child;
        }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::node::get_right_child() const
        {
            return this->right_child;
        }

        template <typename value_t>
        inline bool binary_tree<value_t>::node::has_left_child() const
        {
            return this->left_child != nullptr;
        }

        template <typename value_t>
        inline bool binary_tree<value_t>::node::has_right_child() const
        {
            return this->right_child != nullptr;
        }

        template <typename value_t>
        inline std::vector<std::shared_ptr<typename tree<value_t>::node>> binary_tree<value_t>::node::get_children() const
        {
            std::vector<std::shared_ptr<typename tree<value_t>::node>> children;
            
            if (this->left_child != nullptr)
            {
                children.push_back(std::dynamic_pointer_cast<typename tree<value_t>::node>(this->left_child));
            }

            if (this->right_child != nullptr)
            {
                children.push_back(std::dynamic_pointer_cast<typename tree<value_t>::node>(this->right_child));
            }

            return children;
        }

        template <typename value_t>
        inline std::size_t binary_tree<value_t>::node::get_num_children() const
        {
            return static_cast<std::size_t>(((this->left_child == nullptr) ? 0 : 1) + ((this->right_child == nullptr) ? 0 : 1));
        }

        template <typename value_t>
        inline binary_tree<value_t>::binary_tree() : tree<value_t>(nullptr) { }

        template <typename value_t>
        inline binary_tree<value_t>::binary_tree(std::shared_ptr<node> root) : tree<value_t>(root) { }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::set_root(std::shared_ptr<node> root)
        {
            tree<value_t>::set_root(std::dynamic_pointer_cast<typename tree<value_t>::node>(root));
            return root;
        }

        template <typename value_t>
        inline std::shared_ptr<typename binary_tree<value_t>::node> binary_tree<value_t>::get_root() const
        {
            return std::dynamic_pointer_cast<node>(tree<value_t>::get_root());
        }
    }
}