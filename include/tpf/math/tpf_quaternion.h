#pragma once

#include "Eigen/Dense"

#include <iostream>
#include <type_traits>

namespace tpf
{
    namespace math
    {
        /// <summary>
        /// Quaternion for rotation
        /// </summary>
        /// <template name="floatp_t">Floating point type as value type</template>
        template <typename floatp_t>
        class quaternion
        {
        public:
            using value_type = floatp_t;

            static_assert(std::is_floating_point<floatp_t>::value, "Quaternion must consist of floating point values");

            /// <summary>
            /// Constructor for quaternion representing no rotation (1, [0,0,0])
            /// </summary>
            quaternion();

            /// <summary>
            /// Constructor for quaternion (real, [0,0,0])
            /// </summary>
            /// <param name="real">Real part</param>
            quaternion(floatp_t real);

            /// <summary>
            /// Constructor for quaternion (1, imaginary)
            /// </summary>
            /// <param name="imaginary">Imaginary part</param>
            quaternion(const Eigen::Matrix<floatp_t, 3, 1>& imaginary);

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="real">Real part</param>
            /// <param name="imaginary">Imaginary part</param>
            quaternion(floatp_t real, const Eigen::Matrix<floatp_t, 3, 1>& imaginary);

            /// <summary>
            /// Copy-constructor
            /// </summary>
            /// <param name="copy">Quaternion to copy values from</param>
            quaternion(const quaternion& copy);

            /// <summary>
            /// Assignment operator
            /// </summary>
            /// <param name="copy">Quaternion to copy values from</param>
            /// <returns>This</returns>
            quaternion& operator=(const quaternion& copy);

            /// <summary>
            /// Assignment operator for the real part
            /// </summary>
            /// <param name="new_real">New real part</param>
            /// <returns>This</returns>
            quaternion& operator=(floatp_t new_real);

            /// <summary>
            /// Assignment operator for the imaginary part
            /// </summary>
            /// <param name="new_imaginary">New imaginary part</param>
            /// <returns>This</returns>
            quaternion& operator=(const Eigen::Matrix<floatp_t, 3, 1>& new_imaginary);

            /// <summary>
            /// Multiplicate quaternions
            /// </summary>
            /// <param name="other">Quaternion to multiply with</param>
            /// <returns>Multiplied quaternion</returns>
            quaternion operator*(const quaternion& other) const;

            /// <summary>
            /// Multiplicate quaternions and store result
            /// </summary>
            /// <param name="other">Quaternion to multiply with</param>
            /// <returns>Multiplied quaternion (this)</returns>
            quaternion& operator*=(const quaternion& other);

            /// <summary>
            /// Rotate vector
            /// </summary>
            /// <param name="vec">Vector to rotate</param>
            /// <returns>Rotated vector</returns>
            /// <template name="U">Value type of the vector</template>
            template <typename U>
            Eigen::Matrix<U, 3, 1> operator*(const Eigen::Matrix<U, 3, 1>& vec) const;

            /// <summary>
            /// Equality operator
            /// </summary>
            /// <param name="other">The quaternion to compare with</param>
            /// <returns>True if equal, false otherwise</returns>
            bool operator==(const quaternion& other) const;

            /// <summary>
            /// Get real part
            /// </summary>
            /// <returns>Real part</returns>
            const floatp_t& get_real() const;

            /// <summary>
            /// Get real part
            /// </summary>
            /// <returns>Real part</returns>
            floatp_t& get_real();

            /// <summary>
            /// Set real part
            /// </summary>
            /// <param name="new_real">New real part</returns>
            void set_real(floatp_t new_real);

            /// <summary>
            /// Get imaginary part
            /// </summary>
            /// <returns>Imaginary part</returns>
            const Eigen::Matrix<floatp_t, 3, 1>& get_imaginary() const;

            /// <summary>
            /// Get imaginary part
            /// </summary>
            /// <returns>Imaginary part</returns>
            Eigen::Matrix<floatp_t, 3, 1>& get_imaginary();

            /// <summary>
            /// Set imaginary part
            /// </summary>
            /// <param name="new_imaginary">New imaginary part</returns>
            void set_imaginary(const Eigen::Matrix<floatp_t, 3, 1>& new_imaginary);

            /// <summary>
            /// Return the inverted quaternion
            /// </summary>
            /// <returns>Inverted quaternion</returns>
            quaternion inverse() const;

            /// <summary>
            /// Invert and return the quaternion
            /// </summary>
            /// <returns>Inverted quaternion (this)</returns>
            quaternion& invert();

            /// <summary>
            /// Return the conjugation of the quaternion
            /// </summary>
            /// <returns>Conjugation of the quaternion</returns>
            quaternion conjugation() const;

            /// <summary>
            /// Conjugate and return the quaternion
            /// </summary>
            /// <returns>Conjugated quaternion (this)</returns>
            quaternion& conjugate();

            /// <summary>
            /// Return the normalized of the quaternion
            /// </summary>
            /// <returns>Unit quaternion</returns>
            quaternion normalized() const;

            /// <summary>
            /// Normalize and return the quaternion
            /// </summary>
            /// <returns>Unit quaternion (this)</returns>
            quaternion& normalize();

            /// <summary>
            /// Calculate the norm of the quaternion
            /// </summary>
            /// <returns>Norm of the quaternion</returns>
            floatp_t norm() const;

            /// <summary>
            /// Calculate the squared norm of the quaternion
            /// </summary>
            /// <returns>Squared norm of the quaternion</returns>
            floatp_t squared_norm() const;

            /// <summary>
            /// Create rotation axis
            /// </summary>
            /// <returns>Rotation axis</returns>
            Eigen::Matrix<floatp_t, 3, 1> to_axis() const;

            /// <summary>
            /// Create quaternion from axis representation
            /// </summary>
            /// <param name="axis">Rotation axis</param>
            void from_axis(const Eigen::Matrix<floatp_t, 3, 1>& axis);

        private:
            /// Real part
            floatp_t real;

            /// Imaginary part
            Eigen::Matrix<floatp_t, 3, 1> imaginary;
        };

        /// <summary>
        /// Create quaternion by adding a scalar and a vector
        /// </summary>
        /// <param name="real">Real part</returns>
        /// <param name="imaginary">Imaginary part</returns>
        /// <template name="floatp_t">Value type of the vector and the resulting quaternion</template>
        template <typename floatp_t>
        quaternion<floatp_t> operator+(floatp_t real, const Eigen::Matrix<floatp_t, 3, 1>& imaginary);
        
        /// <summary>
        /// Output quaternion to stream
        /// </summary>
        /// <param name="stream">Output stream</returns>
        /// <param name="quaternion">Quaternion to output</returns>
        /// <template name="floatp_t">Value type of the quaternion</template>
        /// <returns>Output stream</returns>
        template <typename floatp_t>
        std::ostream& operator<<(std::ostream& stream, const quaternion<floatp_t>& quaternion);
    }
}

#include "tpf_quaternion.inl"
