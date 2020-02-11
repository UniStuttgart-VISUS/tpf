#pragma once

#include "Eigen/Dense"

#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Grid information
        /// </summary>
        /// <template name="point_t">Point type</template>
        template <typename point_t>
        struct grid_information
        {
            using value_type = point_t;
            using array_type = std::vector<std::vector<point_t>>;

            array_type cell_coordinates, node_coordinates, cell_sizes;
        };

        /// Convenience for grid coordinates
        using coords1_t = Eigen::Matrix<long long, 1, 1>;
        using coords2_t = Eigen::Matrix<long long, 2, 1>;
        using coords3_t = Eigen::Matrix<long long, 3, 1>;
    }
}