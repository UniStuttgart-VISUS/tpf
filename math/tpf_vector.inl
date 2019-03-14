#include "tpf_vector.h"

#include "Eigen/Dense"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace tpf
{
    namespace math
    {
        template <typename floatp_t>
        inline std::pair<Eigen::Matrix<floatp_t, 3, 1>, Eigen::Matrix<floatp_t, 3, 1>> orthonormal(const Eigen::Matrix<floatp_t, 3, 1>& vector, const floatp_t epsilon)
        {
#ifdef __tpf_sanity_checks
            if (vector.isZero(epsilon))
            {
                throw std::runtime_error("Cannot find orthogonal vectors to a zero-vector");
            }
#endif

            std::pair<Eigen::Matrix<floatp_t, 3, 1>, Eigen::Matrix<floatp_t, 3, 1>> orthonormal_directions;

            // Create first orthogonal direction
            if (vector[0] == static_cast<floatp_t>(0.0L))
            {
                orthonormal_directions.first[0] = static_cast<floatp_t>(0.0L);
                orthonormal_directions.first[1] = -vector[2];
                orthonormal_directions.first[2] = vector[1];
            }
            else if (vector[1] == static_cast<floatp_t>(0.0L))
            {
                orthonormal_directions.first[0] = -vector[2];
                orthonormal_directions.first[1] = static_cast<floatp_t>(0.0L);
                orthonormal_directions.first[2] = vector[0];
            }
            else
            {
                orthonormal_directions.first[0] = -vector[1];
                orthonormal_directions.first[1] = vector[0];
                orthonormal_directions.first[2] = static_cast<floatp_t>(0.0L);
            }

            // Calculate second orthogonal vector using the cross product
            orthonormal_directions.second = vector.cross(orthonormal_directions.first);

            // Normalize and return found vectors
            orthonormal_directions.first.normalize();
            orthonormal_directions.second.normalize();

            return orthonormal_directions;
        }

        template <typename floatp_t>
        inline floatp_t calculate_angle(const Eigen::Matrix<floatp_t, 3, 1>& first, const Eigen::Matrix<floatp_t, 3, 1>& second, const floatp_t epsilon)
        {
#ifdef __tpf_sanity_checks
            if (first.isZero(epsilon) || second.isZero(epsilon))
            {
                throw std::runtime_error("Cannot calculate angle for zero-vectors");
            }
#endif

            return std::acos(first.normalized().dot(second.normalized()));
        }

        template <typename floatp_t, int components>
        inline std::array<Eigen::Matrix<floatp_t, components, 1>, components> get_fitting_directions(const Eigen::Matrix<floatp_t, components, 1>& vector)
        {
            static_assert(components > 0, "Number of vector components must be larger than zero");

            // Sort vector components according to their value
            std::array<std::pair<int, floatp_t>, components> values;

            for (int i = 0; i < components; ++i)
            {
                values[i] = std::make_pair(i, vector[i]);
            }

            std::sort(values.begin(), values.end(),
                [](const std::pair<int, floatp_t>& first, const std::pair<int, floatp_t>& second) { return std::abs(first.second) > std::abs(second.second); });

            // Create directions
            std::array<Eigen::Matrix<floatp_t, components, 1>, components> directions;

            for (int i = 0; i < components; ++i)
            {
                directions[i].setZero();
                directions[i][values[i].first] = std::copysign(static_cast<floatp_t>(1.0), values[i].second);
            }

            return directions;
        }

        template <typename floatp_t, int components = 3>
        inline Eigen::Matrix<floatp_t, components, 1> get_fitting_direction(const Eigen::Matrix<floatp_t, 3, 1>& vector)
        {
            static_assert(components > 0, "Number of vector components must be larger than zero");

            return get_fitting_directions(vector)[0];
        }

        template <typename floatp_t, int components>
        inline bool equals(const Eigen::Matrix<floatp_t, components, 1>& lhs, const Eigen::Matrix<floatp_t, components, 1>& rhs, const floatp_t epsilon)
        {
            bool equal = true;

            for (int i = 0; i < components; ++i)
            {
                equal &= equals(lhs[i], rhs[i], epsilon);
            }

            return equal;
        }
    }
}
