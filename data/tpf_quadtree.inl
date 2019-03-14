#include "tpf_quadtree.h"

#include "tpf_n-tree.h"
#include "tpf_position.h"
#include "tpf_tree.h"

#include "../geometry/tpf_containing.h"
#include "../geometry/tpf_point.h"
#include "../geometry/tpf_rectangle.h"

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
        inline quadtree<float_t, value_t, kernel_t>::node::node(const store_type& value) : base_type::node(value)
        { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline quadtree<float_t, value_t, kernel_t>::node::node(store_type&& value) : base_type::node(std::move(value))
        { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline auto quadtree<float_t, value_t, kernel_t>::node::set_child(const position_t position, const value_t& value) -> std::shared_ptr<typename base_type::node>
        {
            // Split area of current node and calculate new node's part
            const auto rectangle = tree<store_type>::node::get_value().first;

            const auto diagonal = rectangle->get_max_point().get_internal() - rectangle->get_min_point().get_internal();
            const auto center = rectangle->get_min_point().get_internal() + 0.5 * diagonal;

            const typename kernel_t::Vector_3 offset_x(0.5 * diagonal.x(), 0.0, 0.0);
            const typename kernel_t::Vector_3 offset_y(0.0, 0.5 * diagonal.y(), 0.0);

            std::shared_ptr<geometry::rectangle<float_t, kernel_t>> new_area = nullptr;

            switch (position)
            {
            case position_t::BOTTOM_LEFT:
                new_area = std::make_shared<geometry::rectangle<float_t, kernel_t>>(rectangle->get_min_point().get_internal(), center);
                break;
            case position_t::BOTTOM_RIGHT:
                new_area = std::make_shared<geometry::rectangle<float_t, kernel_t>>(rectangle->get_min_point().get_internal() + offset_x, center + offset_x);
                break;
            case position_t::TOP_LEFT:
                new_area = std::make_shared<geometry::rectangle<float_t, kernel_t>>(rectangle->get_min_point().get_internal() + offset_y, center + offset_y);
                break;
            case position_t::TOP_RIGHT:
                new_area = std::make_shared<geometry::rectangle<float_t, kernel_t>>(center, rectangle->get_max_point().get_internal());
            }

            // Create new node
            auto child = std::make_shared<node>(std::make_pair(new_area, value));

            child->set_parent(this);

            this->children[static_cast<std::size_t>(position)] = child;

            return child;
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline void quadtree<float_t, value_t, kernel_t>::node::set_children(const std::array<value_t, 4>& values)
        {
            set_child(position_t::BOTTOM_LEFT, values[0]);
            set_child(position_t::BOTTOM_RIGHT, values[1]);
            set_child(position_t::TOP_LEFT, values[2]);
            set_child(position_t::TOP_RIGHT, values[3]);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline quadtree<float_t, value_t, kernel_t>::quadtree() { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline quadtree<float_t, value_t, kernel_t>::quadtree(std::shared_ptr<node> root) : base_type(std::dynamic_pointer_cast<typename base_type::node>(root)) { }

        template <typename float_t, typename value_t, typename kernel_t>
        inline auto quadtree<float_t, value_t, kernel_t>::find_node(const geometry::point<float_t, kernel_t>& point) const ->
            std::pair<std::shared_ptr<typename base_type::node>, std::list<position_t>>
        {
            // Find leaf node in which the point is located
            auto current = base_type::get_root();

            if (!geometry::is_in_rectangle(point, *current->get_value().first))
            {
                throw std::runtime_error(__tpf_error_message("Point is located outside of the quadtree"));
            }

            std::list<position_t> path;

            while (true)
            {
                // Get children
                auto bottom_left = current->get_child(position_t::BOTTOM_LEFT);
                auto bottom_right = current->get_child(position_t::BOTTOM_RIGHT);
                auto top_left = current->get_child(position_t::TOP_LEFT);
                auto top_right = current->get_child(position_t::TOP_RIGHT);

                // Locate point in child nodes and add turn taken to path
                if (bottom_left != nullptr && geometry::is_in_rectangle(point, *bottom_left->get_value().first))
                {
                    current = bottom_left;
                    path.push_back(position_t::BOTTOM_LEFT);
                }
                else if (bottom_right != nullptr && geometry::is_in_rectangle(point, *bottom_right->get_value().first))
                {
                    current = bottom_right;
                    path.push_back(position_t::BOTTOM_RIGHT);
                }
                else if (top_left != nullptr && geometry::is_in_rectangle(point, *top_left->get_value().first))
                {
                    current = top_left;
                    path.push_back(position_t::TOP_LEFT);
                }
                else if (top_right != nullptr && geometry::is_in_rectangle(point, *top_right->get_value().first))
                {
                    current = top_right;
                    path.push_back(position_t::TOP_RIGHT);
                }
                else
                {
                    return std::make_pair(current, path);
                }
            }
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline value_t quadtree<float_t, value_t, kernel_t>::interpolate_at_center(const std::vector<std::shared_ptr<typename tree<store_type>::node>>& nodes) const
        {
            return math::interpolate_bilinear(static_cast<float_t>(0.5), static_cast<float_t>(0.5),
                nodes[0]->get_value().second, nodes[1]->get_value().second, nodes[2]->get_value().second, nodes[3]->get_value().second);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline value_t quadtree<float_t, value_t, kernel_t>::interpolate_at_point(const geometry::point<float_t, kernel_type>& point,
            const std::vector<std::pair<geometry::point<float_t, kernel_t>, value_t>>& nodes) const
        {
            return math::interpolate_bilinear(point.get_vertex(), nodes[0].first.get_vertex(), nodes[1].first.get_vertex(),
                nodes[2].first.get_vertex(), nodes[3].first.get_vertex(), nodes[0].second, nodes[1].second, nodes[2].second, nodes[3].second);
        }

        template <typename float_t, typename value_t, typename kernel_t>
        inline std::array<neighbor_t, 3> quadtree<float_t, value_t, kernel_t>::get_directions(const geometry::point<float_t, kernel_type>& relative_point) const
        {
            std::array<neighbor_t, 3> directions = {
                relative_point.get_vertex()[0] > 0.0 ? neighbor_t::RIGHT : neighbor_t::LEFT,
                relative_point.get_vertex()[1] > 0.0 ? neighbor_t::TOP : neighbor_t::BOTTOM,
                neighbor_t::NONE
            };

            directions[2] = static_cast<neighbor_t>(static_cast<int>(directions[0]) | static_cast<int>(directions[1]));

            return directions;
        }
    }
}
