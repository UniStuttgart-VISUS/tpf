#include "tpf_octree.h"

#include "../data/tpf_octree.h"
#include "../data/tpf_position.h"

#include "../geometry/tpf_cuboid.h"
#include "../geometry/tpf_point.h"

#include "../log/tpf_log.h"

#include "vtkDataArray.h"

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        template <typename point_t, typename value_t, int num_components>
        inline data::octree<point_t, octree_aux::value_type<value_t, num_components>>
            get_octree(vtkDataArray* id, std::function<octree_aux::value_type<value_t, num_components>(vtkIdType)> get_values, const std::array<double, 6>& bounds)
        {
            using stored_t = octree_aux::value_type<value_t, num_components>;
            using octree_t = data::octree<point_t, stored_t>;

            if (id->GetNumberOfComponents() != 1)
            {
                throw std::runtime_error(__tpf_error_message("Number of components must be 1 for octree node paths."));
            }

            // Create octree with root using the bounding box from the data
            const geometry::point<point_t> min_point(bounds[0], bounds[2], bounds[4]);
            const geometry::point<point_t> max_point(bounds[1], bounds[3], bounds[5]);

            auto bounding_box = std::make_shared<geometry::cuboid<point_t>>(min_point, max_point);
            auto root = std::make_shared<typename octree_t::node>(std::make_pair(bounding_box, stored_t()));

            octree_t octree(root);

            // Add nodes
            std::vector<std::pair<std::vector<data::position_t>, stored_t>> node_information;
            node_information.reserve(id->GetNumberOfValues());

            for (vtkIdType p = 0; p < id->GetNumberOfValues(); ++p)
            {
                // Decode path
                auto path_code = id->GetVariantValue(p).ToLongLong();

                std::vector<data::position_t> path(static_cast<std::size_t>(std::floor(std::log2(path_code) / std::log2(8))), data::position_t::INVALID);
                std::size_t index = 0;

                while (path_code != 1)
                {
                    const auto turn = path_code % 8;

                    path[index++] = static_cast<data::position_t>(((turn & 1) ? 4 : 0) + ((turn & 2) ? 2 : 0) + ((turn & 4) ? 1 : 0));

                    path_code /= 8;
                }

                // Add node information
                node_information.push_back(std::make_pair(path, get_values(p)));
            }

            octree.insert_nodes(node_information.begin(), node_information.end());

            return octree;
        }
    }
}
