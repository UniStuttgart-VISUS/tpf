#pragma once

#include "../math/tpf_transformer.h"

#include "Eigen/Dense"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <memory>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// Geometry types
        enum class geometry_t
        {
            POINT, LINE, TRIANGLE, TETRAHEDRON, PLANE, RECTANGLE, CUBOID, POLYGON, POLYHEDRON
        };

        /// Default CGAL kernel
        using default_kernel_t = CGAL::Exact_predicates_inexact_constructions_kernel;

        /// <summary>
        /// Base class for geometric objects
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        template <typename floatp_t>
        class geometric_object
        {
        public:
            static_assert(std::is_floating_point<floatp_t>::value, "Must be a floating point type");

            using value_type = floatp_t;

            /// <summary>
            /// Clone object (deep copy)
            /// </summary>
            /// <returns>Deep copy</returns>
            virtual std::shared_ptr<geometric_object> clone(const math::transformer<floatp_t, 3>& trafo) const = 0;

            /// <summary>
            /// Transform this object using a transformer
            /// </summary>
            /// <param name="trafo">Transformer</param>
            virtual geometric_object<floatp_t>& transform(const math::transformer<floatp_t, 3>& trafo) = 0;

            /// <summary>
            /// Return number of points
            /// </summary>
            /// <returns>Number of points</returns>
            virtual std::size_t get_num_points() const = 0;

            /// <summary>
            /// Return points
            /// </summary>
            /// <returns>Points</returns>
            virtual std::vector<Eigen::Matrix<floatp_t, 3, 1>> get_points() const = 0;

            /// <summary>
            /// Return number of cells
            /// </summary>
            /// <returns>Number of cells</returns>
            virtual std::size_t get_num_cells() const = 0;

            /// <summary>
            /// Return cells
            /// </summary>
            /// <returns>Cells</returns>
            virtual std::vector<std::vector<std::size_t>> get_cells() const = 0;

            /// <summary>
            /// Serialize object to a binary representation
            /// </summary>
            /// <returns>Binary representation</returns>
            virtual std::vector<char> serialize() const = 0;
        };
    }
}
