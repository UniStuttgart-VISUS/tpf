#pragma once

#include "tpf_quaternion.h"

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
            /// Transformation to a coordinate system defined by an up-point and a normal
            /// </summary>
            /// <param name="point">Origin</param>
            /// <param name="normal">Normal</param>
            /// <param name="invert">Invert transformation matrix</param>
            transformer(const Eigen::Matrix<floatp_t, 3, 1>& point, const Eigen::Matrix<floatp_t, 3, 1>& normal, bool invert = false);

            /// <summary>
            /// Transformation to a given coordinate system
            /// </summary>
            /// <param name="origin">Origin</param>
            /// <param name="x_axis">x axis</param>
            /// <param name="y_axis">y axis</param>
            /// <param name="z_axis">z axis</param>
            /// <param name="invert">Invert transformation matrix</param>
            transformer(const Eigen::Matrix<floatp_t, 3, 1>& origin, const Eigen::Matrix<floatp_t, 3, 1>& x_axis,
                const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis, bool invert = false);

            /// <summary>
            /// Transformation from translation and rotation
            /// </summary>
            /// <param name="translation">Translation</param>
            /// <param name="quaternion">Rotation quaternion</param>
            /// <param name="center_of_rotation">Center around to rotate</param>
            /// <param name="invert">Invert transformation matrix</param>
            transformer(const Eigen::Matrix<floatp_t, 3, 1>& translation, const quaternion<floatp_t>& quaternion,
                const Eigen::Matrix<floatp_t, 3, 1>& center_of_rotation, bool invert = false);

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
            /// Multiply internal transformation matrices
            /// </summary>
            /// <param name="other">Other transformer</param>
            /// <returns>Transformer containing the product of the transformation matrices</returns>
            transformer operator*(const transformer& other) const;

            /// <summary>
            /// Multiply internal transformation matrices
            /// </summary>
            /// <param name="other">Other transformer</param>
            /// <returns>This, containing the product of the transformation matrices</returns>
            transformer& operator*=(const transformer& other);

            /// <summary>
            /// Answer if this transformer would not alter the input
            /// </summary>
            /// <returns>True if the transformer does not alter the input; false otherwise</returns>
            bool is_unit() const;

            /// <summary>
            /// Set pre-processing function called before transformation
            /// </summary>
            /// <param name="func">Pre-processing function</param>
            void set_preprocessing(std::function<Eigen::Matrix<floatp_t, 4, 1>(const Eigen::Matrix<floatp_t, 4, 1>&)> func);

            /// <summary>
            /// Invert this transformation
            /// </summary>
            /// <returns>This, containing the inverse transformation</returns>
            transformer& invert();

            /// <summary>
            /// Inverse transformation
            /// </summary>
            /// <returns>Inverted transformation</returns>
            transformer inverse() const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform(const Eigen::Matrix<floatp_t, dimension, 1>& vec, floatp_t w = homogeneous) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform_inverse(const Eigen::Matrix<floatp_t, dimension, 1>& vec, floatp_t w = homogeneous) const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec, floatp_t w = homogeneous) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_inverse_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec, floatp_t w = homogeneous) const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform_normal(const Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1> transform_normal_inverse(const Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_normal_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

            /// <summary>
            /// Transform vector inverse
            /// </summary>
            /// <param name="vec">Vector to transform</param>
            /// <return>Transformed vector</return>
            Eigen::Matrix<floatp_t, dimension, 1>& transform_normal_inverse_inplace(Eigen::Matrix<floatp_t, dimension, 1>& vec) const;

        private:
            /// <summary>
            /// Create matrix for transformation from input into standard coordinate system
            /// </summary>
            /// <param name="origin">Origin of input coordinate system</param>
            /// <param name="x_axis">X axis of input coordinate system</param>
            /// <param name="y_axis">Y axis of input coordinate system</param>
            /// <param name="z_axis">Z axis of input coordinate system</param>
            /// <param name="invert">Invert transformation matrix</param>
            void create_transformation_matrix(const Eigen::Matrix<floatp_t, 3, 1>& origin, const Eigen::Matrix<floatp_t, 3, 1>& x_axis,
                const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis, bool invert);

            /// <summary>
            /// Set default pre-processing function
            /// </summary>
            void set_default_preprocessing() noexcept;

            /// Transformation matrix
            Eigen::Matrix<floatp_t, 4, 4> trafo;

            /// Pre-processing function
            std::function<Eigen::Matrix<floatp_t, 4, 1>(const Eigen::Matrix<floatp_t, 4, 1>&)> preprocessing;
            bool is_default_function;
        };
    }
}

#include "tpf_transformer.inl"
