#include "tpf_octree.h"

#include "tpf_n-tree.h"
#include "tpf_position.h"
#include "tpf_tree.h"

#include "../geometry/tpf_containing.h"
#include "../geometry/tpf_cuboid.h"
#include "../geometry/tpf_point.h"

#include "../math/tpf_interpolation.h"

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
        template <typename float_t, typename value_t, typename kernel_t>
        inline octree<float_t, value_t, kernel_t>::node::node(const store_type& value) : base_type::node(value)
        { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline octree<float_t, value_t, kernel_t>::node::node(store_type&& value) : base_type::node(std::move(value))
        { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline auto octree<float_t, value_t, kernel_t>::node::set_child(const position_t position, const value_t& value) -> std::shared_ptr<typename base_type::node>
        {
            // Split area of current node and calculate new node's part
            const auto cuboid = tree<store_type>::node::get_value().first;

            const auto diagonal = cuboid->get_max_point().get_internal() - cuboid->get_min_point().get_internal();
            const auto center = cuboid->get_min_point().get_internal() + 0.5 * diagonal;

            const typename kernel_t::Vector_3 offset_x(0.5 * diagonal.x(), 0.0, 0.0);
            const typename kernel_t::Vector_3 offset_y(0.0, 0.5 * diagonal.y(), 0.0);
            const typename kernel_t::Vector_3 offset_z(0.0, 0.0, 0.5 * diagonal.z());

            std::shared_ptr<geometry::cuboid<float_t, kernel_t>> new_volume = nullptr;

            switch (position)
            {
            case position_t::FRONT_BOTTOM_LEFT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal(), center);
                break;
            case position_t::FRONT_BOTTOM_RIGHT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_x, center + offset_x);
                break;
            case position_t::FRONT_TOP_LEFT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_y, center + offset_y);
                break;
            case position_t::FRONT_TOP_RIGHT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_x + offset_y, center + offset_x + offset_y);
                break;
            case position_t::BACK_BOTTOM_LEFT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_z, center + offset_z);
                break;
            case position_t::BACK_BOTTOM_RIGHT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_x + offset_z, center + offset_x + offset_z);
                break;
            case position_t::BACK_TOP_LEFT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(cuboid->get_min_point().get_internal() + offset_y + offset_z, center + offset_y + offset_z);
                break;
            case position_t::BACK_TOP_RIGHT:
                new_volume = std::make_shared<geometry::cuboid<float_t, kernel_t>>(center, cuboid->get_max_point().get_internal());
            }

            // Create new node
            auto child = std::make_shared<node>(std::make_pair(new_volume, value));

            child->set_parent(this);

            this->children[static_cast<std::size_t>(position)] = child;

            return child;
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline void octree<float_t, value_t, kernel_t>::node::set_children(const std::array<value_t, 8>& values)
        {
            set_child(position_t::FRONT_BOTTOM_LEFT, values[0]);
            set_child(position_t::FRONT_BOTTOM_RIGHT, values[1]);
            set_child(position_t::FRONT_TOP_LEFT, values[2]);
            set_child(position_t::FRONT_TOP_RIGHT, values[3]);
            set_child(position_t::BACK_BOTTOM_LEFT, values[4]);
            set_child(position_t::BACK_BOTTOM_RIGHT, values[5]);
            set_child(position_t::BACK_TOP_LEFT, values[6]);
            set_child(position_t::BACK_TOP_RIGHT, values[7]);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline octree<float_t, value_t, kernel_t>::octree() { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline octree<float_t, value_t, kernel_t>::octree(std::shared_ptr<node> root) : base_type(std::dynamic_pointer_cast<typename base_type::node>(root)) { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline auto octree<float_t, value_t, kernel_t>::find_node(const geometry::point<float_t, kernel_t>& point) const ->
            std::pair<std::shared_ptr<typename base_type::node>, std::list<position_t>>
        {
            // Find leaf node in which the point is located
            auto current = base_type::get_root();

            if (!geometry::is_in_cuboid(point, *current->get_value().first))
            {
                throw std::runtime_error(__tpf_error_message("Point is located outside of the octree"));
            }

            std::list<position_t> path;

            while (true)
            {
                // Get children
                auto front_bottom_left = current->get_child(position_t::FRONT_BOTTOM_LEFT);
                auto front_bottom_right = current->get_child(position_t::FRONT_BOTTOM_RIGHT);
                auto front_top_left = current->get_child(position_t::FRONT_TOP_LEFT);
                auto front_top_right = current->get_child(position_t::FRONT_TOP_RIGHT);
                auto back_bottom_left = current->get_child(position_t::BACK_BOTTOM_LEFT);
                auto back_bottom_right = current->get_child(position_t::BACK_BOTTOM_RIGHT);
                auto back_top_left = current->get_child(position_t::BACK_TOP_LEFT);
                auto back_top_right = current->get_child(position_t::BACK_TOP_RIGHT);

                // Locate point in child nodes and add turn taken to path
                if (front_bottom_left != nullptr && geometry::is_in_cuboid(point, *front_bottom_left->get_value().first))
                {
                    current = front_bottom_left;
                    path.push_back(position_t::FRONT_BOTTOM_LEFT);
                }
                else if (front_bottom_right != nullptr && geometry::is_in_cuboid(point, *front_bottom_right->get_value().first))
                {
                    current = front_bottom_right;
                    path.push_back(position_t::FRONT_BOTTOM_RIGHT);
                }
                else if (front_top_left != nullptr && geometry::is_in_cuboid(point, *front_top_left->get_value().first))
                {
                    current = front_top_left;
                    path.push_back(position_t::FRONT_TOP_LEFT);
                }
                else if (front_top_right != nullptr && geometry::is_in_cuboid(point, *front_top_right->get_value().first))
                {
                    current = front_top_right;
                    path.push_back(position_t::FRONT_TOP_RIGHT);
                }
                else if (back_bottom_left != nullptr && geometry::is_in_cuboid(point, *back_bottom_left->get_value().first))
                {
                    current = back_bottom_left;
                    path.push_back(position_t::BACK_BOTTOM_LEFT);
                }
                else if (back_bottom_right != nullptr && geometry::is_in_cuboid(point, *back_bottom_right->get_value().first))
                {
                    current = back_bottom_right;
                    path.push_back(position_t::BACK_BOTTOM_RIGHT);
                }
                else if (back_top_left != nullptr && geometry::is_in_cuboid(point, *back_top_left->get_value().first))
                {
                    current = back_top_left;
                    path.push_back(position_t::BACK_TOP_LEFT);
                }
                else if (back_top_right != nullptr && geometry::is_in_cuboid(point, *back_top_right->get_value().first))
                {
                    current = back_top_right;
                    path.push_back(position_t::BACK_TOP_RIGHT);
                }
                else
                {
                    return std::make_pair(current, path);
                }
            }
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline value_t octree<float_t, value_t, kernel_t>::interpolate_at_center(const std::vector<std::shared_ptr<typename tree<store_type>::node>>& nodes) const
        {
            return math::interpolate_trilinear(static_cast<float_t>(0.5), static_cast<float_t>(0.5), static_cast<float_t>(0.5),
                nodes[0]->get_value().second, nodes[1]->get_value().second, nodes[2]->get_value().second, nodes[3]->get_value().second,
                nodes[4]->get_value().second, nodes[5]->get_value().second, nodes[6]->get_value().second, nodes[7]->get_value().second);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline value_t octree<float_t, value_t, kernel_t>::interpolate_at_point(const geometry::point<float_t, kernel_type>& point,
            const std::vector<std::pair<geometry::point<float_t, kernel_t>, value_t>>& nodes) const
        {
            return math::interpolate_trilinear(point.get_vertex(), nodes[0].first.get_vertex(), nodes[1].first.get_vertex(),
                nodes[2].first.get_vertex(), nodes[3].first.get_vertex(), nodes[4].first.get_vertex(), nodes[5].first.get_vertex(),
                nodes[6].first.get_vertex(), nodes[7].first.get_vertex(), nodes[0].second, nodes[1].second, nodes[2].second,
                nodes[3].second, nodes[4].second, nodes[5].second, nodes[6].second, nodes[7].second);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline std::array<neighbor_t, 7> octree<float_t, value_t, kernel_t>::get_directions(const geometry::point<float_t, kernel_type>& relative_point) const
        {
            std::array<neighbor_t, 7> directions = {
                relative_point.get_internal()[0] > 0.0 ? neighbor_t::RIGHT : neighbor_t::LEFT,
                relative_point.get_internal()[1] > 0.0 ? neighbor_t::TOP : neighbor_t::BOTTOM,
                neighbor_t::NONE,
                relative_point.get_internal()[2] > 0.0 ? neighbor_t::BACK : neighbor_t::FRONT,
                neighbor_t::NONE,
                neighbor_t::NONE,
                neighbor_t::NONE
            };

            directions[2] = static_cast<neighbor_t>(static_cast<int>(directions[0]) | static_cast<int>(directions[1]));
            directions[4] = static_cast<neighbor_t>(static_cast<int>(directions[0]) | static_cast<int>(directions[3]));
            directions[5] = static_cast<neighbor_t>(static_cast<int>(directions[1]) | static_cast<int>(directions[3]));
            directions[6] = static_cast<neighbor_t>(static_cast<int>(directions[0]) | static_cast<int>(directions[1]) | static_cast<int>(directions[3]));

            return directions;
        }
    }
}
