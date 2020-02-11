#include "tpf_linalg.h"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

#include <cmath>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace math
    {
        template <typename value_t, typename point_t, std::size_t dimension>
        inline Eigen::Matrix<value_t, dimension, 1> calculate_gradient(const std::vector<Eigen::Matrix<point_t, dimension, 1>>& vertices,
            const std::vector<value_t>& values, const bool inverse_distance_weighting)
        {
            static_assert(std::is_floating_point<value_t>::value, "Value type must be floating point");
            static_assert(std::is_floating_point<point_t>::value, "Point type must be floating point");

            Eigen::Matrix<value_t, dimension, 1> gradient;

            if (vertices.size() < 3 || vertices.size() != values.size())
            {
                gradient.setZero();
                return gradient;
            }

            // Create matrix A, weight matrix W and right-hand side b
            Eigen::Matrix<value_t, Eigen::Dynamic, dimension> A;
            A.resize(vertices.size() - 1, Eigen::NoChange);

            Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> W;
            W.resize(vertices.size() - 1, vertices.size() - 1);
            W.setZero();

            Eigen::Matrix<value_t, Eigen::Dynamic, 1> b;
            b.resize(vertices.size() - 1, Eigen::NoChange);

            for (std::size_t j = 0; j < vertices.size() - 1; ++j)
            {
                for (std::size_t d = 0; d < dimension; ++d)
                {
                    A(j, d) = static_cast<value_t>(vertices[j + 1][d] - vertices[0][d]);
                }

                if (inverse_distance_weighting)
                {
                    W(j, j) = static_cast<value_t>(1.0L) / (vertices[j + 1] - vertices[0]).squaredNorm();
                }
                else
                {
                    W(j, j) = static_cast<value_t>(1.0L);
                }

                b(j) = values[j + 1] - values[0];
            }

            // Calculate weighted matrix and vectors
            Eigen::Matrix<value_t, Eigen::Dynamic, dimension> X;
            X.resize(vertices.size() - 1, Eigen::NoChange);

            X = W * A;
            b = W * b;

            // Solve least squares for Ax = b
            Eigen::JacobiSVD<Eigen::Matrix<value_t, Eigen::Dynamic, dimension>> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);

            gradient = svd.solve(b);

            return gradient;
        }

        template <typename value_t, typename point_t, std::size_t dimension>
        inline Eigen::Matrix<value_t, dimension, dimension> calculate_jacobian(const std::vector<Eigen::Matrix<point_t, dimension, 1>>& vertices,
            const std::vector<Eigen::Matrix<value_t, dimension, 1>>& values, const bool inverse_distance_weighting)
        {
            static_assert(std::is_floating_point<value_t>::value, "Value type must be floating point");
            static_assert(std::is_floating_point<point_t>::value, "Point type must be floating point");

            Eigen::Matrix<value_t, dimension, dimension> jacobian;

            if (vertices.size() < 3 || vertices.size() != values.size())
            {
                jacobian.setZero();
                return jacobian;
            }

            // Get values component-wise
            auto x = std::vector<value_t>(values.size());
            auto y = std::vector<value_t>(values.size());
            auto z = std::vector<value_t>(values.size());

            for (std::size_t i = 0; i < values.size(); ++i)
            {
                x[i] = values[i][0];
                y[i] = values[i][1];
                z[i] = values[i][2];
            }

            // Solve least squares for Ax = b
            Eigen::Matrix<value_t, dimension, 1> result[3];

            result[0] = calculate_gradient<value_t, point_t, dimension>(vertices, x, inverse_distance_weighting);
            result[1] = calculate_gradient<value_t, point_t, dimension>(vertices, y, inverse_distance_weighting);
            result[2] = calculate_gradient<value_t, point_t, dimension>(vertices, z, inverse_distance_weighting);

            // Return jacobian
            for (std::size_t i = 0; i < 3; i++)
            {
                for (std::size_t j = 0; j < 3; j++)
                {
                    jacobian(j, i) = result[j][i];
                }
            }

            return jacobian;
        }

        template <typename floatp_t, std::size_t rows_and_columns>
        inline Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> extract_symmetric_part(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix)
        {
            static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

            Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> symmetric = static_cast<floatp_t>(0.5L) * (matrix + matrix.transpose());

            return symmetric;
        }

        template <typename floatp_t, std::size_t rows_and_columns>
        inline Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> extract_asymmetric_part(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix)
        {
            static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

            Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns> asymmetric = static_cast<floatp_t>(0.5L) * (matrix - matrix.transpose());

            return asymmetric;
        }

        namespace
        {
            template <typename floatp_t, std::size_t rows_and_columns>
            inline eigenpairs<floatp_t, rows_and_columns> calculate_eigenpair_symm(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix)
            {
                static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

                eigenpairs<floatp_t, rows_and_columns> ret;

                const floatp_t a = matrix(0, 0);
                const floatp_t b = matrix(0, 1);
                const floatp_t c = matrix(1, 1);

                ret.eigenvalues[0] = (a + c + std::sqrt((a - c) * (a - c) + static_cast<floatp_t>(4.0) * b * b)) / static_cast<floatp_t>(2.0);
                ret.eigenvalues[1] = (a + c - std::sqrt((a - c) * (a - c) + static_cast<floatp_t>(4.0) * b * b)) / static_cast<floatp_t>(2.0);

                ret.eigenvectors[0][0] = static_cast<floatp_t>(1.0);
                ret.eigenvectors[0][1] = (ret.eigenvalues[0] - a) / b;
                ret.eigenvectors[0].normalize();

                ret.eigenvectors[1][0] = ret.eigenvectors[0][1];
                ret.eigenvectors[1][1] = -ret.eigenvectors[0][0];

                return ret;
            }

            template <typename floatp_t, std::size_t rows_and_columns>
            inline eigenpairs<floatp_t, rows_and_columns> calculate_eigenpair_asymm(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix)
            {
                static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

                Eigen::EigenSolver<Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>> eigensolver(matrix, true);

                eigenpairs<floatp_t, rows_and_columns> ret;

                for (std::size_t index = 0; index < rows_and_columns; ++index)
                {
                    ret.eigenvectors[index] = eigensolver.eigenvectors().real().col(index);
                    ret.eigenvalues[index] = eigensolver.eigenvalues().real()[index];
                }

                return ret;
            }
        }

        template <typename floatp_t, std::size_t rows_and_columns>
        inline eigenpairs<floatp_t, rows_and_columns> calculate_eigenpair(const Eigen::Matrix<floatp_t, rows_and_columns, rows_and_columns>& matrix)
        {
            if (rows_and_columns == 2 && matrix(0, 1) == matrix(1, 0))
            {
                return calculate_eigenpair_symm<floatp_t, rows_and_columns>(matrix);
            }
            else
            {
                return calculate_eigenpair_asymm<floatp_t, rows_and_columns>(matrix);
            }
        }
    }
}