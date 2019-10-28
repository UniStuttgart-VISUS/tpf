#pragma once

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace math
    {
        /// <summary>
        /// Transformation of vectors
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="dimension">Vector dimension</template>
        /// <template name="homogeneous">homogeneous coordinate for 3d vectors</template>
        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous = 1>
        class transformer
        {
        public:
            static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

            /// <summary>
            /// Get transformer storing a unit matrix
            /// </summary>
            /// <returns>Unit transformer</returns>
            static transformer unit();

            /// <summary>
            /// Constructor
            /// </summary>
            transformer();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="point">Origin</param>
            /// <param name="normal">Normal</param>
            transformer(const Eigen::Matrix<floatp_t, 3, 1>& point, const Eigen::Matrix<floatp_t, 3, 1>& normal);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="origin">Origin</param>
            /// <param name="x_axis">x axis</param>
            /// <param name="y_axis">y axis</param>
            /// <param name="z_axis">z axis</param>
            transformer(const Eigen::Matrix<floatp_t, 3, 1>& origin, const Eigen::Matrix<floatp_t, 3, 1>& x_axis,
                const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="trafo">Transformation matrix</param>
            transformer(const Eigen::Matrix<floatp_t, 4, 4>& trafo);

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Transformer to copy from</param>
            transformer(const transformer& copy);

            /// <summary>
            /// Move constructor
            /// </summary>
            /// <param name="move">Transformer to move from</param>
            transformer(transformer&& move);

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Transformer to copy from</param>
            transformer& operator=(const transformer& copy);

            /// <summary>
            /// Move assignment
            /// </summary>
            /// <param name="move">Transformer to move from</param>
            transformer& operator=(transformer&& move);

            /// <summary>
            /// Check equality of transformations
            /// </summary>
            /// <param name="rhs">Other transformer to compare with</param>
            /// <returns>True if same transformer; false otherwise</returns>
            bool operator==(const transformer& rhs) const;

            /// <summary>
            /// Check inequality of transformations
            /// </summary>
            /// <param name="rhs">Other transformer to compare with</param>
            /// <returns>True if not the same transformer; false otherwise</returns>
            bool operator!=(const transformer& rhs) const;

            /// <summary>
            /// Answer if this transformer stores a unit matrix
            /// </summary>
            /// <returns>Does this transformer store a unit matrix?</returns>
            bool is_unit() const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform(const Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform_inverse(const Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_inverse_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

        private:
            /// <summary>
            /// Create matrix for transformation from input into standard coordinate system
            /// </summary>
            /// <param name="origin">Origin of input coordinate system</param>
            /// <param name="x_axis">X axis of input coordinate system</param>
            /// <param name="y_axis">Y axis of input coordinate system</param>
            /// <param name="z_axis">Z axis of input coordinate system</param>
            void create_transformation_matrix(const Eigen::Matrix<floatp_t, 3, 1>& origin, const Eigen::Matrix<floatp_t, 3, 1>& x_axis,
                const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis);

            /// Transformation matrix
            Eigen::Matrix<floatp_t, 4, 4> trafo;
        };
    }
}

#include "tpf_transformer.inl"