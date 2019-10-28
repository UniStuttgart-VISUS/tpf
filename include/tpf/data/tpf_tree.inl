#include "tpf_tree.h"

#include <list>
#include <memory>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename value_t>
        inline tree<value_t>::node::node()
        {
            this->parent = nullptr;
        }

        template <typename value_t>
        inline tree<value_t>::node::~node()
        { }

        template <typename value_t>
        inline bool tree<value_t>::node::is_leaf() const
        {
            return get_num_children() == 0;
        }

        template <typename value_t>
        inline const value_t& tree<value_t>::node::get_value() const
        {
            return this->value;
        }

        template <typename value_t>
        inline value_t& tree<value_t>::node::get_value()
        {
            return this->value;
        }

        template <typename value_t>
        inline void tree<value_t>::node::set_parent(node* parent)
        {
            this->parent = parent;
        }

        template <typename value_t>
        inline typename tree<value_t>::node* tree<value_t>::node::get_parent() const
        {
            return this->parent;
        }

        template <typename value_t>
        inline tree<value_t>::tree()
        {
            this->root = nullptr;
        }

        template <typename value_t>
        inline tree<value_t>::tree(std::shared_ptr<node> root)
        {
            this->root = root;
        }

        template <typename value_t>
        inline std::shared_ptr<typename tree<value_t>::node> tree<value_t>::set_root(std::shared_ptr<node> root)
        {
            this->root = root;
            return root;
        }

        template <typename value_t>
        inline std::shared_ptr<typename tree<value_t>::node> tree<value_t>::get_root() const
        {
            return this->root;
        }

        template <typename value_t>
        inline std::vector<std::list<std::shared_ptr<typename tree<value_t>::node>>> tree<value_t>::get_paths() const
        {
            return get_paths(this->root);
        }

        template <typename value_t>
        inline std::vector<std::list<std::shared_ptr<typename tree<value_t>::node>>> tree<value_t>::get_paths(std::shared_ptr<node> root) const
        {
            std::vector<std::list<std::shared_ptr<node>>> paths;

            if (root->is_leaf())
            {
                paths.push_back(std::list<std::shared_ptr<node>> { root });
                
                return paths;
            }

            for (std::size_t i = 0; i < root->get_num_children(); ++i)
            {
                auto subpaths = get_paths(root->get_children()[i]);

                for (auto& subpath : subpaths)
                {
                    subpath.push_front(root);
                }

                paths.insert(paths.end(), subpaths.begin(), subpaths.end());
            }

            return paths;
        }

        template <typename value_t>
        inline std::list<std::shared_ptr<typename tree<value_t>::node>> tree<value_t>::get_leaf_nodes() const
        {
            auto nodes = depth_first_search(this->root);

            for (auto it = nodes.begin(); it != nodes.end();)
            {
                if (!(*it)->is_leaf())
                {
                    nodes.erase(it++);
                }
                else
                {
                    ++it;
                }
            }

            return nodes;
        }

        template <typename value_t>
        inline std::list<std::shared_ptr<typename tree<value_t>::node>> tree<value_t>::breadth_first_search() const
        {
            std::list<std::shared_ptr<node>> queue;
            queue.push_back(this->root);

            std::list<std::shared_ptr<node>> nodes;

            while (!queue.empty())
            {
                auto node = queue.front();
                queue.pop_front();

                nodes.push_back(node);

                auto children = node->get_children();

                queue.insert(queue.end(), children.begin(), children.end());
            }

            return nodes;
        }

        template <typename value_t>
        inline std::list<std::shared_ptr<typename tree<value_t>::node>> tree<value_t>::depth_first_search() const
        {
            return depth_first_search(this->root);
        }

        template <typename value_t>
        inline std::list<std::shared_ptr<typename tree<value_t>::node>> tree<value_t>::depth_first_search(std::shared_ptr<node> root) const
        {
            std::list<std::shared_ptr<typename tree<value_t>::node>> nodes;

            nodes.push_back(root);

            if (root->is_leaf())
            {
                return nodes;
            }

            for (std::size_t i = 0; i < root->get_num_children(); ++i)
            {
                auto nodes_from_subtrees = depth_first_search(root->get_children()[i]);

                nodes.insert(nodes.end(), nodes_from_subtrees.begin(), nodes_from_subtrees.end());
            }

            return nodes;
        }
    }
}