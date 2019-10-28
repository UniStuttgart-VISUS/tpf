#include "tpf_bsp_tree.h"

#include "tpf_binary_tree.h"
#include "tpf_tree.h"

#include "../geometry/tpf_halfspace.h"
#include "../geometry/tpf_plane.h"
#include "../geometry/tpf_point.h"

#include "../math/tpf_math_defines.h"
#include "../math/tpf_vector.h"

#include "Eigen/Dense"

#include <list>
#include <memory>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename float_t, typename kernel_t>
        inline bsp_tree<float_t, kernel_t>::bsp_tree() : binary_tree<value_type>(nullptr) { }

        template <typename float_t, typename kernel_t>
        inline bsp_tree<float_t, kernel_t>::bsp_tree(std::shared_ptr<node_type> root) : binary_tree<value_type>(root) { }

        template <typename float_t, typename kernel_t>
        template <typename forward_iterator_t>
        inline void bsp_tree<float_t, kernel_t>::create_bsp_tree(forward_iterator_t begin, forward_iterator_t end)
        {
            // Create root node
            auto root_partition = std::make_shared<geometry::plane<float_t, kernel_t>>(construct_partitioning_plane(begin->first->get_points(), begin->second));

            auto node = std::make_shared<node_type>(std::make_pair(root_partition, begin->second));
            binary_tree<value_type>::set_root(node);

            // Initialize values for the root node
            std::list<std::pair<std::shared_ptr<node_type>, std::list<input_value_type>>> node_and_values;

            node_and_values.push_back(std::make_pair(node, std::list<input_value_type>(++begin, end)));

            // For each node, split it according to the partitioning plane
            while (!node_and_values.empty())
            {
                const auto& key_value = node_and_values.front();

                const auto node = key_value.first;
                const auto& list = key_value.second;

                // Create partitioning plane and sort the values into two sets, based on their relative position
                std::list<input_value_type> in, out;

                const auto partition = node->get_value().first;

                for (auto it = list.begin(); it != list.end(); ++it)
                {
                    const auto input_points = it->first->get_points();

                    bool inside, outside;
                    inside = outside = false;

                    for (const auto& point : input_points)
                    {
                        inside |= geometry::is_in_negative_halfspace(geometry::point<float_t, kernel_t>(point), *partition);
                        outside |= geometry::is_in_positive_halfspace(geometry::point<float_t, kernel_t>(point), *partition);
                    }

                    if (inside)
                    {
                        in.push_back(*it);
                    }
                    if (outside)
                    {
                        out.push_back(*it);
                    }
                }

                // Add nodes and their values to the set of nodes to be handled
                if (in.size() > 0)
                {
                    auto in_partition = construct_partitioning_plane(in.front().first->get_points(), in.front().second);
                    auto in_node = std::make_shared<node_type>(std::make_pair(std::make_shared<geometry::plane<float_t, kernel_t>>(in_partition), in.front().second));
                    in.pop_front();

                    node->set_left_child(in_node);

                    if (in.size() > 0)
                    {
                        node_and_values.push_back(std::make_pair(in_node, in));
                    }
                }

                if (out.size() > 0)
                {
                    auto out_partition = construct_partitioning_plane(out.front().first->get_points(), out.front().second);
                    auto out_node = std::make_shared<node_type>(std::make_pair(std::make_shared<geometry::plane<float_t, kernel_t>>(out_partition), out.front().second));
                    out.pop_front();

                    node->set_right_child(out_node);

                    if (out.size() > 0)
                    {
                        node_and_values.push_back(std::make_pair(out_node, out));
                    }
                }

                node_and_values.pop_front();
            }
        }

        template <typename float_t, typename kernel_t>
        inline bool bsp_tree<float_t, kernel_t>::is_inside(const geometry::point<float_t, kernel_t>& point) const
        {
            std::shared_ptr<node_type> node = std::dynamic_pointer_cast<node_type>(tree<value_type>::get_root());

            while (true)
            {
                // Get next node or, if possible, decide directly if inside or not
                if (geometry::is_on_plane(point, *node->get_value().first))
                {
                    return true;
                }

                if (geometry::is_in_negative_halfspace(point, *node->get_value().first))
                {
                    if (node->has_left_child())
                    {
                        node = node->get_left_child();
                    }
                    else
                    {
                        return true;
                    }
                }
                else
                {
                    if (node->has_right_child())
                    {
                        node = node->get_right_child();
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }

        template <typename float_t, typename kernel_t>
        inline std::vector<bsp_tree<float_t, kernel_t>> bsp_tree<float_t, kernel_t>::extract_convex_subtrees() const
        {
            const auto paths = tree<value_type>::get_paths();
            
            std::vector<bsp_tree<float_t, kernel_t>> subtrees;
            subtrees.reserve(paths.size());

            // For each path, create a new tree
            for (const auto& path : paths)
            {
                std::shared_ptr<node_type> last_node = nullptr;
                bool last_node_was_left_child = true;

                for (auto it = path.rbegin(); it != path.rend(); ++it)
                {
                    auto node_copy = std::make_shared<node_type>((*it)->get_value());

                    if (last_node != nullptr)
                    {
                        if (last_node_was_left_child)
                        {
                            node_copy->set_left_child(last_node);
                        }
                        else
                        {
                            node_copy->set_right_child(last_node);
                        }
                    }

                    last_node = node_copy;

                    if ((*it)->get_parent() != nullptr)
                    {
                        if (dynamic_cast<node_type*>((*it)->get_parent())->get_left_child() == *it)
                        {
                            last_node_was_left_child = true;
                        }
                        else
                        {
                            last_node_was_left_child = false;
                        }
                    }
                }

                subtrees.emplace_back(last_node);
            }

            return subtrees;
        }

        template <typename float_t, typename kernel_t>
        inline geometry::plane<float_t, kernel_t> bsp_tree<float_t, kernel_t>::construct_partitioning_plane(
            const std::vector<Eigen::Matrix<float_t, 3, 1>>& points, const Eigen::Matrix<float_t, 3, 1>& normal) const
        {
            const geometry::point<float_t, kernel_t> first(points[0]);
            const geometry::point<float_t, kernel_t> second(points[1]);
            const geometry::point<float_t, kernel_t> third(points[2]);

            const Eigen::Matrix<float_t, 3, 1> calculated_normal = (points[1] - points[0]).cross(points[2] - points[0]).normalized();

            if (calculated_normal.isApprox(normal) || math::calculate_angle(calculated_normal, normal) < (math::pi<float_t> / static_cast<float_t>(2.0)))
            {
                return geometry::plane<float_t, kernel_t>(first, second, third);
            }
            else
            {
                return geometry::plane<float_t, kernel_t>(first, third, second);
            }
        }
    }
}