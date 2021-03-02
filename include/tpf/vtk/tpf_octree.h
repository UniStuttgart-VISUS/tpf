#pragma once

#include "../data/tpf_octree.h"

#include "vtkDataArray.h"

#include <array>
#include <functional>
#include <type_traits>

namespace tpf
{
    namespace vtk
    {
        namespace octree_aux
        {
            template <typename value_t, int num_components>
            using value_type = typename std::conditional<num_components == 1, value_t, Eigen::Matrix<value_t, num_components, 1>>::type;
        }

        /// <summary>
        /// Construct an octree from node IDs and stored values
        /// </summary>
        /// <typeparam name="point_t">Point type</typeparam>
        /// <typeparam name="value_t">Type of stored value</typeparam>
        /// <template name="num_components">Number of components of the value vector</template>
        /// <param name="id">ID indicating the path within the octree</param>
        /// <param name="get_values">Get values stored at the nodes</param>
        /// <param name="bounds">Boundaries of the spatial domain</param>
        /// <returns>Constructed octree</returns>
        template <typename point_t, typename value_t, int num_components>
        data::octree<point_t, octree_aux::value_type<value_t, num_components>> get_octree(vtkDataArray* id,
            std::function<octree_aux::value_type<value_t, num_components>(vtkIdType)> get_values, const std::array<double, 6>& bounds);
    }
}

#include "tpf_octree.inl"
