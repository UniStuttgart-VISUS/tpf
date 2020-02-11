#include "tpf_n-tree.h"

#include "tpf_position.h"
#include "tpf_tree.h"

#include "../exception/tpf_not_implemented_exception.h"

#include "../geometry/tpf_point.h"

#include "../log/tpf_log.h"

#include "Eigen/Dense"

#include <algorithm>
#include <array>
#include <list>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        namespace
        {
            template <typename number_t>
            struct initialize_zero
            {
                inline number_t operator()()
                {
                    return static_cast<number_t>(0);
                }
            };

            template <typename number_t, int rows, int columns>
            struct initialize_zero<Eigen::Matrix<number_t, rows, columns>>
            {
                inline Eigen::Matrix<number_t, rows, columns> operator()()
                {
                    Eigen::Matrix<number_t, rows, columns> vector;
                    vector.setZero();

                    return vector;
                }
            };
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline n_tree<float_t, value_t, geometry_t, num_children>::node::node(const store_type& value)
        {
            tree<store_type>::node::get_value() = value;

            this->children.fill(nullptr);
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline n_tree<float_t, value_t, geometry_t, num_children>::node::node(store_type&& value)
        {
            tree<store_type>::node::get_value() = std::move(value);

            this->children.fill(nullptr);
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline std::shared_ptr<typename n_tree<float_t, value_t, geometry_t, num_children>::node> n_tree<float_t, value_t, geometry_t, num_children>::node::get_child(const position_t position) const
        {
            return this->children[static_cast<std::size_t>(position)];
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline bool n_tree<float_t, value_t, geometry_t, num_children>::node::has_child(const position_t position) const
        {
            return this->children[static_cast<std::size_t>(position)] != nullptr;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline auto n_tree<float_t, value_t, geometry_t, num_children>::node::get_children() const -> std::vector<std::shared_ptr<typename tree<store_type>::node>>
        {
            std::vector<std::shared_ptr<typename tree<store_type>::node>> children;

            for (const auto child : this->children)
            {
                if (child != nullptr)
                {
                    children.push_back(std::dynamic_pointer_cast<typename tree<store_type>::node>(child));
                }
            }

            return children;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline std::size_t n_tree<float_t, value_t, geometry_t, num_children>::node::get_num_children() const
        {
            return std::count_if(this->children.begin(), this->children.end(), [](const std::shared_ptr<node> child) { return child != nullptr; });
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline n_tree<float_t, value_t, geometry_t, num_children>::n_tree() { }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline n_tree<float_t, value_t, geometry_t, num_children>::n_tree(std::shared_ptr<node> root) : tree<store_type>(std::dynamic_pointer_cast<typename tree<store_type>::node>(root)) { }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline std::shared_ptr<typename n_tree<float_t, value_t, geometry_t, num_children>::node> n_tree<float_t, value_t, geometry_t, num_children>::set_root(std::shared_ptr<node> root)
        {
            tree<store_type>::set_root(std::dynamic_pointer_cast<typename tree<store_type>::node>(root));
            return root;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline std::shared_ptr<typename n_tree<float_t, value_t, geometry_t, num_children>::node> n_tree<float_t, value_t, geometry_t, num_children>::get_root() const
        {
            return std::dynamic_pointer_cast<node>(tree<store_type>::get_root());
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        template <typename iterator_t>
        void n_tree<float_t, value_t, geometry_t, num_children>::insert_nodes(iterator_t begin, iterator_t end)
        {
            // Insert all nodes
            for (iterator_t current = begin; current != end; ++current)
            {
                insert_node(current->first, current->second);
            }

            // Check for well-formedness of tree, i.e. that every node has zero or four children
#ifdef __tpf_sanity_checks
            const auto nodes = tree<store_type>::depth_first_search();

            for (const auto node : nodes)
            {
                if (node->get_num_children() != 0 && node->get_num_children() != num_children)
                {
                    throw std::runtime_error(__tpf_error_message("Tree is not well-formed. All nodes of an n-tree must have either zero or n children."));
                }
            }
#endif
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline std::shared_ptr<typename n_tree<float_t, value_t, geometry_t, num_children>::node> n_tree<float_t, value_t, geometry_t, num_children>
            ::insert_node(const std::vector<position_t>& position, value_t& value)
        {
            // Find node to which the new one should be attached to
            auto last_node = get_root();
            auto current_node = last_node;

            std::size_t index = 0;

            for (; index < position.size() - 1 && current_node != nullptr; ++index)
            {
                last_node = current_node;
                current_node = last_node->get_child(position[index]);
            }

            // Create intermediate nodes until the given position can be inserted
            if (current_node == nullptr)
            {
                --index;
                current_node = last_node;
            }

            // Create new node(s) with given value
            for (; index < position.size(); ++index)
            {
                current_node = current_node->set_child(position[index], value);
            }

            // If inserted node was the last missing child, interpolate value for all parents containing four nodes
            auto current_parent = current_node.get();

            while (current_parent->get_parent() != nullptr && current_parent->get_parent()->get_num_children() == num_children)
            {
                current_parent = dynamic_cast<node*>(current_parent->get_parent());

                current_parent->get_value().second = interpolate_at_center(current_parent->get_children());
            }

            return current_node;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        std::list<position_t> n_tree<float_t, value_t, geometry_t, num_children>::get_node_path(const std::shared_ptr<node> _node) const
        {
            std::list<position_t> path;

            auto* child_node = _node.get();
            auto* current_node = _node.get();

            while (current_node->get_parent() != nullptr)
            {
                current_node = static_cast<node*>(current_node->get_parent());

                for (unsigned int i = 0; i < num_children; ++i)
                {
                    if (current_node->get_children()[i].get() == child_node)
                    {
                        path.push_front(static_cast<position_t>(i));
                        break;
                    }
                }

                child_node = current_node;
            }

            return path;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline auto n_tree<float_t, value_t, geometry_t, num_children>::get_neighbors(const std::shared_ptr<node> center) const -> std::vector<std::shared_ptr<node>>
        {
            const auto directions = get_directions();

            // Compute hypothetical paths for each direction
            const auto node_path = get_node_path(center);
            std::array<std::vector<position_t>, num_children - 1> neighbor_paths;

            for (std::size_t neighbor_index = 0; neighbor_index < num_children - 1; ++neighbor_index)
            {
                neighbor_paths[neighbor_index].insert(neighbor_paths[neighbor_index].end(), node_path.begin(), node_path.end());

                // Find deepest possible turn
                for (auto it = neighbor_paths[neighbor_index].rbegin(); it != neighbor_paths[neighbor_index].rend(); ++it)
                {
                    const auto possible_position = *it + directions[neighbor_index];

                    if (possible_position != position_t::INVALID)
                    {
                        *it = possible_position;

                        // Adjust all deeper positions
                        if (it != neighbor_paths[neighbor_index].rbegin())
                        {
                            for (--it; it != neighbor_paths[neighbor_index].rbegin(); --it)
                            {
                                *it = *it - directions[neighbor_index];
                            }

                            *it = *it - directions[neighbor_index];
                        }

                        break;
                    }
                }
            }

            // Check for existance of neighboring nodes
            std::vector<std::shared_ptr<node>> neighbors;
            neighbors.reserve(num_children);

            for (auto& neighbor_path : neighbor_paths)
            {
                auto neighbor_node = get_node(neighbor_path);

                // Move higher up the tree until the requested node exists or the root has been reached
                while (neighbor_node == nullptr && !neighbor_path.empty())
                {
                    neighbor_path.pop_back();

                    neighbor_node = get_node(neighbor_path);
                }

                // If there is a valid neighbor node, store it
                if (neighbor_node != nullptr)
                {
                    neighbors.push_back(neighbor_node);
                }
            }

            return neighbors;
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline value_t n_tree<float_t, value_t, geometry_t, num_children>::interpolate(const geometry::point<float_t, kernel_type>& point) const
        {
            try
            {
                // Find node in which the point is located
                const auto deepest_point_node = find_node(point);

                // Starting at the deepest level...
                std::array<std::shared_ptr<node>, num_children - 1> neighbor_nodes;
                auto used_point_node = std::make_pair<node*, std::list<position_t>>(deepest_point_node.first.get(), std::list<position_t>(deepest_point_node.second));

                bool success = false;

                while (used_point_node.first != nullptr && !success)
                {
                    // ...get directions
                    const auto point_cell = used_point_node.first->get_value().first;

                    const auto mid_point = 0.5 * ((point_cell->get_max_point().get_internal() - CGAL::ORIGIN) + (point_cell->get_min_point().get_internal() - CGAL::ORIGIN));
                    const auto relative_point = point.get_internal() - mid_point;

                    const auto directions = get_directions(geometry::point<float_t, kernel_type>(static_cast<float_t>(relative_point.x()),
                        static_cast<float_t>(relative_point.y()), static_cast<float_t>(relative_point.z())));

                    // Compute path for each direction
                    std::array<std::vector<position_t>, num_children - 1> neighbor_paths;

                    success = true;

                    for (std::size_t neighbor_index = 0; neighbor_index < num_children - 1 && success; ++neighbor_index)
                    {
                        neighbor_paths[neighbor_index].insert(neighbor_paths[neighbor_index].end(), used_point_node.second.begin(), used_point_node.second.end());

                        // Find deepest possible turn
                        for (auto it = neighbor_paths[neighbor_index].rbegin(); it != neighbor_paths[neighbor_index].rend(); ++it)
                        {
                            const auto possible_position = *it + directions[neighbor_index];

                            if (possible_position != position_t::INVALID)
                            {
                                *it = possible_position;

                                // Adjust all deeper positions
                                if (it != neighbor_paths[neighbor_index].rbegin())
                                {
                                    for (--it; it != neighbor_paths[neighbor_index].rbegin(); --it)
                                    {
                                        *it = *it - directions[neighbor_index];
                                    }

                                    *it = *it - directions[neighbor_index];
                                }

                                break;
                            }
                        }

                        success = !std::equal(neighbor_paths[neighbor_index].begin(), neighbor_paths[neighbor_index].end(), used_point_node.second.begin());
                    }

                    // Check for existance of neighbor nodes
                    if (success)
                    {
                        for (std::size_t neighbor_index = 0; neighbor_index < num_children - 1; ++neighbor_index)
                        {
                            neighbor_nodes[neighbor_index] = get_node(neighbor_paths[neighbor_index]);

                            if (neighbor_nodes[neighbor_index] == nullptr)
                            {
                                // Move higher up the tree and try again
                                used_point_node.first = dynamic_cast<node*>(used_point_node.first->get_parent());

                                if (used_point_node.second.size() > 0)
                                {
                                    used_point_node.second.pop_back();
                                }

                                success = false;

                                break;
                            }
                        }
                    }
                    else
                    {
                        // Move higher up the tree and try again
                        used_point_node.first = dynamic_cast<node*>(used_point_node.first->get_parent());

                        if (used_point_node.second.size() > 0)
                        {
                            used_point_node.second.pop_back();
                        }

                        success = false;
                    }
                }

                if (!success)
                {
#ifdef __tpf_debug
                    log::warning_message(__tpf_warning_message("Interpolation failed and returned 0-value. Point is probably located outside of the interpolatable region"));
#endif

                    return initialize_zero<value_t>()();
                }

                // Get points and respective values for interpolation
                std::vector<std::pair<geometry::point<float_t, kernel_type>, value_t>> points;

                {
                    const auto point_cell = used_point_node.first->get_value().first;
                    const auto mid_point = 0.5 * ((point_cell->get_max_point().get_internal() - CGAL::ORIGIN) + (point_cell->get_min_point().get_internal() - CGAL::ORIGIN));

                    points.push_back(std::make_pair(geometry::point<float_t, kernel_type>(static_cast<float_t>(mid_point.x()),
                        static_cast<float_t>(mid_point.y()), static_cast<float_t>(mid_point.z())), used_point_node.first->get_value().second));
                }

                for (const auto& neighbor_node : neighbor_nodes)
                {
                    const auto point_cell = neighbor_node->get_value().first;
                    const auto mid_point = 0.5 * ((point_cell->get_max_point().get_internal() - CGAL::ORIGIN) + (point_cell->get_min_point().get_internal() - CGAL::ORIGIN));

                    points.push_back(std::make_pair(geometry::point<float_t, kernel_type>(static_cast<float_t>(mid_point.x()),
                        static_cast<float_t>(mid_point.y()), static_cast<float_t>(mid_point.z())), neighbor_node->get_value().second));
                }

                // Interpolate the value at the given position
                return interpolate_at_point(point, points);
            }
            catch (const std::runtime_error&
#ifdef __tpf_debug
                e
#endif
                )
            {
#ifdef __tpf_debug
                log::warning_message(__tpf_nested_warning_message(e.what(), "Interpolation failed and returned 0-value. Point is probably located outside of the covered domain"));
#endif

                return initialize_zero<value_t>()();
            }

#ifdef __tpf_debug
            log::warning_message(__tpf_warning_message("Interpolation failed and returned 0-value"));
#endif

            return initialize_zero<value_t>()();
        }

        template <typename float_t, typename value_t, typename geometry_t, std::size_t num_children>
        inline auto n_tree<float_t, value_t, geometry_t, num_children>::get_node(const std::vector<position_t>& position) const -> std::shared_ptr<node>
        {
            std::shared_ptr<node> requested_node = std::dynamic_pointer_cast<node>(tree<store_type>::get_root());

            for (auto pos : position)
            {
                if (requested_node == nullptr || pos == position_t::INVALID)
                {
                    return nullptr;
                }

                requested_node = requested_node->get_child(pos);
            }

            return requested_node;
        }
    }
}
