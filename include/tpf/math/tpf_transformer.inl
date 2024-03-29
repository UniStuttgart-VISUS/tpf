#include "tpf_transformer.h"

#include "tpf_vector.h"

#include "Eigen/Dense"

#include <utility>

namespace tpf
{
    namespace math
    {
        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous> transformer<floatp_t, dimension, homogeneous>::unit()
        {
            return transformer<floatp_t, dimension, homogeneous>();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer()
        {
            this->trafo.setIdentity();
            set_default_preprocessing();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(const Eigen::Matrix<floatp_t, 3, 1>& point,
            const Eigen::Matrix<floatp_t, 3, 1>& normal, const bool invert)
        {
            const auto directions = orthonormal<floatp_t>(normal);

            create_transformation_matrix(point, directions.first, directions.second, normal, invert);
            set_default_preprocessing();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(const Eigen::Matrix<floatp_t, 3, 1>& origin, const Eigen::Matrix<floatp_t, 3, 1>& x_axis,
            const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis, const bool invert)
        {
            create_transformation_matrix(origin, x_axis, y_axis, z_axis, invert);
            set_default_preprocessing();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(const Eigen::Matrix<floatp_t, 3, 1>& translation,
            const quaternion<floatp_t>& quaternion, const Eigen::Matrix<floatp_t, 3, 1>& center_of_rotation, const bool invert)
        {
            Eigen::Matrix<floatp_t, 4, 4> translate;
            translate.setIdentity();
            translate.col(3).head(3) = translation;

            Eigen::Matrix<floatp_t, 4, 4> to_origin;
            to_origin.setIdentity();
            to_origin.col(3).head(3) = -center_of_rotation;

            Eigen::Matrix<floatp_t, 4, 4> to_center;
            to_center.setIdentity();
            to_center.col(3).head(3) = center_of_rotation;

            Eigen::Matrix<floatp_t, 4, 4> rotate;
            rotate.setIdentity();
            rotate.block(0, 0, 3, 3) = quaternion.to_matrix();

            this->trafo = translate * to_center * rotate * to_origin;
            set_default_preprocessing();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(const Eigen::Matrix<floatp_t, 4, 4>& trafo)
        {
            this->trafo = trafo;
            set_default_preprocessing();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(const transformer& copy)
        {
            this->trafo = copy.trafo;
            this->preprocessing = copy.preprocessing;
            this->is_default_function = copy.is_default_function;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>::transformer(transformer&& move)
        {
            this->trafo = std::move(move.trafo);
            this->preprocessing = move.preprocessing;
            this->is_default_function = move.is_default_function;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>& transformer<floatp_t, dimension, homogeneous>::operator=(const transformer& copy)
        {
            this->trafo = copy.trafo;
            this->preprocessing = copy.preprocessing;
            this->is_default_function = copy.is_default_function;

            return *this;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>& transformer<floatp_t, dimension, homogeneous>::operator=(transformer&& move)
        {
            this->trafo = std::move(move.trafo);
            this->preprocessing = move.preprocessing;
            this->is_default_function = move.is_default_function;

            return *this;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline bool transformer<floatp_t, dimension, homogeneous>::operator==(const transformer& rhs) const
        {
            return this->trafo == rhs.trafo;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous> transformer<floatp_t, dimension, homogeneous>::operator*(const transformer& other) const
        {
            return transformer<floatp_t, dimension, homogeneous>(this->trafo * other.trafo);
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>& transformer<floatp_t, dimension, homogeneous>::operator*=(const transformer& other)
        {
            this->trafo *= other.trafo;

            return *this;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline bool transformer<floatp_t, dimension, homogeneous>::is_unit() const
        {
            return (this->trafo == unit().trafo) && this->is_default_function;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline void transformer<floatp_t, dimension, homogeneous>::set_preprocessing(std::function<Eigen::Matrix<floatp_t, 4, 1>(
            const Eigen::Matrix<floatp_t, 4, 1>&)> func)
        {
            this->preprocessing = func;
            this->is_default_function = false;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous>& transformer<floatp_t, dimension, homogeneous>::invert()
        {
            this->trafo = this->trafo.inverse();
            return *this;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline transformer<floatp_t, dimension, homogeneous> transformer<floatp_t, dimension, homogeneous>::inverse() const
        {
            transformer copy(*this);
            return copy.invert();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1> transformer<floatp_t, dimension, homogeneous>::transform(
            const Eigen::Matrix<floatp_t, dimension, 1>& vec, const floatp_t w) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, w;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo * this->preprocessing(temp);

            // Return
            return temp.template head<dimension>();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1> transformer<floatp_t, dimension, homogeneous>::transform_inverse(
            const Eigen::Matrix<floatp_t, dimension, 1>& vec, const floatp_t w) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, w;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.inverse() * this->preprocessing(temp);

            // Return
            return temp.template head<dimension>();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1>& transformer<floatp_t, dimension, homogeneous>::transform_inplace(
            Eigen::Matrix<floatp_t, dimension, 1>& vec, const floatp_t w) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, w;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo * this->preprocessing(temp);

            // Return
            vec << temp.template head<dimension>();

            return vec;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1>& transformer<floatp_t, dimension, homogeneous>::transform_inverse_inplace(
            Eigen::Matrix<floatp_t, dimension, 1>& vec, const floatp_t w) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, w;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.inverse() * this->preprocessing(temp);

            // Return
            vec << temp.template head<dimension>();

            return vec;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1> transformer<floatp_t, dimension, homogeneous>::transform_normal(
            const Eigen::Matrix<floatp_t, dimension, 1>& vec) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, 0.0;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.inverse().transpose() * temp;

            // Return
            return temp.template head<dimension>();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1> transformer<floatp_t, dimension, homogeneous>::transform_normal_inverse(
            const Eigen::Matrix<floatp_t, dimension, 1>& vec) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, 0.0;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.transpose() * temp;

            // Return
            return temp.template head<dimension>();
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1>& transformer<floatp_t, dimension, homogeneous>::transform_normal_inplace(
            Eigen::Matrix<floatp_t, dimension, 1>& vec) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, 0.0;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.inverse().transpose() * temp;

            // Return
            vec << temp.template head<dimension>();

            return vec;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline Eigen::Matrix<floatp_t, dimension, 1>& transformer<floatp_t, dimension, homogeneous>::transform_normal_inverse_inplace(
            Eigen::Matrix<floatp_t, dimension, 1>& vec) const
        {
            static_assert(dimension == 3 || dimension == 4, "Dimension must be either 3 or 4");

            // Calculate transformation
            Eigen::Matrix<floatp_t, 4, 1> temp;

            if (dimension == 3)
            {
                temp << vec, 0.0;
            }
            else if (dimension == 4)
            {
                temp << vec;
            }

            temp = this->trafo.transpose() * temp;

            // Return
            vec << temp.template head<dimension>();

            return vec;
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline void transformer<floatp_t, dimension, homogeneous>::create_transformation_matrix(const Eigen::Matrix<floatp_t, 3, 1>& origin,
            const Eigen::Matrix<floatp_t, 3, 1>& x_axis, const Eigen::Matrix<floatp_t, 3, 1>& y_axis, const Eigen::Matrix<floatp_t, 3, 1>& z_axis,
            const bool invert)
        {
            this->trafo <<
                x_axis[0], y_axis[0], z_axis[0], origin[0],
                x_axis[1], y_axis[1], z_axis[1], origin[1],
                x_axis[2], y_axis[2], z_axis[2], origin[2],
                static_cast<floatp_t>(0.0L), static_cast<floatp_t>(0.0L), static_cast<floatp_t>(0.0L), static_cast<floatp_t>(1.0L);

            if (invert)
            {
                this->trafo = this->trafo.inverse();
            }
        }

        template <typename floatp_t, std::size_t dimension, std::size_t homogeneous>
        inline void transformer<floatp_t, dimension, homogeneous>::set_default_preprocessing() noexcept
        {
            this->preprocessing = [](const Eigen::Matrix<floatp_t, 4, 1>& vector) { return vector; };
            this->is_default_function = true;
        }
    }
}
