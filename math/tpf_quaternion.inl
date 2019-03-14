#include "tpf_quaternion.h"

#include "tpf_math_defines.h"

#include "Eigen/Dense"

#include <cmath>
#include <iostream>

namespace tpf
{
    namespace math
    {
        template <typename floatp_t>
        inline quaternion<floatp_t>::quaternion() : real(1.0)
        { }

        template <typename floatp_t>
        inline quaternion<floatp_t>::quaternion(floatp_t real) : real(real)
        { }

        template <typename floatp_t>
        inline quaternion<floatp_t>::quaternion(const Eigen::Matrix<floatp_t, 3, 1>& imaginary) : real(1.0), imaginary(imaginary)
        { }

        template <typename floatp_t>
        inline quaternion<floatp_t>::quaternion(floatp_t real, const Eigen::Matrix<floatp_t, 3, 1>& imaginary) : real(real), imaginary(imaginary)
        { }

        template <typename floatp_t>
        inline quaternion<floatp_t>::quaternion(const quaternion<floatp_t>& copy)
        {
            this->real = copy.real;
            this->imaginary = copy.imaginary;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::operator=(const quaternion<floatp_t>& copy)
        {
            this->real = copy.real;
            this->imaginary = copy.imaginary;

            return *this;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::operator=(floatp_t new_real)
        {
            this->real = new_real;
            this->imaginary.setZero();

            return *this;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::operator=(const Eigen::Matrix<floatp_t, 3, 1>& new_imaginary)
        {
            this->real = 1.0;
            this->imaginary = new_imaginary;

            return *this;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t> quaternion<floatp_t>::operator*(const quaternion<floatp_t>& other) const
        {
            quaternion<floatp_t> copy(*this);
            return copy *= other;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::operator*=(const quaternion<floatp_t>& other)
        {
            const floatp_t real = this->real * other.real - this->imaginary.dot(other.imaginary);
            const Eigen::Matrix<floatp_t, 3, 1> imaginary = this->real * other.imaginary + other.real * this->imaginary + this->imaginary.cross(other.imaginary);

            this->real = real;
            this->imaginary = imaginary;

            return *this;
        }

        template <typename floatp_t>
        template <typename U>
        inline Eigen::Matrix<U, 3, 1> quaternion<floatp_t>::operator*(const Eigen::Matrix<U, 3, 1>& vec) const
        {
            return ((*this) * quaternion<floatp_t>(static_cast<floatp_t>(0.0), vec) * this->inverse()).get_imaginary();
        }

        template <typename floatp_t>
        inline bool quaternion<floatp_t>::operator==(const quaternion& other) const
        {
            return equals(this->real, other.real, static_cast<floatp_t>(std::min(this->real, other.real) * 0.00001)) && this->imaginary.isApprox(other.imaginary);
        }

        template <typename floatp_t>
        inline const floatp_t& quaternion<floatp_t>::get_real() const
        {
            return this->real;
        }

        template <typename floatp_t>
        inline floatp_t& quaternion<floatp_t>::get_real()
        {
            return this->real;
        }

        template <typename floatp_t>
        inline void quaternion<floatp_t>::set_real(floatp_t new_real)
        {
            this->real = new_real;
        }

        template <typename floatp_t>
        inline const Eigen::Matrix<floatp_t, 3, 1>& quaternion<floatp_t>::get_imaginary() const
        {
            return this->imaginary;
        }

        template <typename floatp_t>
        inline Eigen::Matrix<floatp_t, 3, 1>& quaternion<floatp_t>::get_imaginary()
        {
            return this->imaginary;
        }

        template <typename floatp_t>
        inline void quaternion<floatp_t>::set_imaginary(const Eigen::Matrix<floatp_t, 3, 1>& new_imaginary)
        {
            this->imaginary = new_imaginary;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t> quaternion<floatp_t>::inverse() const
        {
            quaternion<floatp_t> copy(*this);
            return copy.invert();
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::invert()
        {
            const floatp_t sqr_norm = squared_norm();

            conjugate();
            this->real /= sqr_norm;
            this->imaginary /= sqr_norm;

            return *this;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t> quaternion<floatp_t>::conjugation() const
        {
            quaternion<floatp_t> copy(*this);
            return copy.conjugate();
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::conjugate()
        {
            this->imaginary *= -1.0;

            return *this;
        }

        template <typename floatp_t>
        inline quaternion<floatp_t> quaternion<floatp_t>::normalized() const
        {
            quaternion<floatp_t> copy(*this);
            return copy.normalize();
        }

        template <typename floatp_t>
        inline quaternion<floatp_t>& quaternion<floatp_t>::normalize()
        {
            const floatp_t my_norm = norm();

            this->real /= my_norm;
            this->imaginary /= my_norm;

            return *this;
        }

        template <typename floatp_t>
        inline floatp_t quaternion<floatp_t>::norm() const
        {
            return std::sqrt(squared_norm());
        }

        template <typename floatp_t>
        inline floatp_t quaternion<floatp_t>::squared_norm() const
        {
            return this->real * this->real + this->imaginary.squaredNorm();
        }

        template <typename floatp_t>
        inline Eigen::Matrix<floatp_t, 3, 1> quaternion<floatp_t>::to_axis() const
        {
            return Eigen::Matrix<floatp_t, 3, 1>(static_cast<floatp_t>(2.0L) * std::atan2(this->imaginary.norm(), this->real) * this->imaginary.normalized());
        }

        template <typename floatp_t>
        inline void quaternion<floatp_t>::from_axis(const Eigen::Matrix<floatp_t, 3, 1>& axis)
        {
            this->real = std::cos(axis.norm() / static_cast<floatp_t>(2.0));
            this->imaginary = axis.normalized() * std::sin(axis.norm() / static_cast<floatp_t>(2.0));
        }

        template <typename floatp_t>
        inline quaternion<floatp_t> operator+(floatp_t real, const Eigen::Matrix<floatp_t, 3, 1>& imaginary)
        {
            return quaternion<floatp_t>(real, imaginary);
        }

        template <typename floatp_t>
        inline std::ostream& operator<<(std::ostream& stream, const quaternion<floatp_t>& quaternion)
        {
            stream << "[" << quaternion.get_real() << ", (" << quaternion.get_imaginary()[0]
                << ", " << quaternion.get_imaginary()[1] << ", " << quaternion.get_imaginary()[2] << ")]";

            return stream;
        }
    }
}
