#pragma once

#include "Eigen/Dense"

#include <vector>

namespace tpf
{
    namespace algorithm
    {
        /// <summary>
        /// Second-order polynomial for use in least squares, storing a polynomial of the form 'a_xx x^2 + a_yy y^2 + a_xy xy + a_x x + a_y y + a_1'
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        template <typename floatp_t>
        struct polynomial
        {
            static_assert(std::is_floating_point<floatp_t>::value, "Type must be floating point");

            /// Stored values for the polynomial of the form 'a_xx x^2 + a_yy y^2 + a_xy xy + a_x x + a_y y + a_1'
            floatp_t a_xx;
            floatp_t a_yy;
            floatp_t a_xy;
            floatp_t a_x;
            floatp_t a_y;
            floatp_t a_1;

            /// <summary>
            /// Subtraction
            /// </summary>
            /// <param name="poly">Other polynomial</param>
            /// <return>Difference</return>
            polynomial<floatp_t> operator-(const polynomial<floatp_t>& poly) const;
        };

        /// <summary>
        /// Least squares for a resulting second-order polynomial
        /// </summary>
        /// <remarks>Needs at least 6 positions, throws an exception when given less</remarks>
        /// <template name="forward_iterator_t">Forward iterator type</template>
        /// <param name="begin">Iterator pointing at the first position vector of type Eigen::Vector3</param>
        /// <param name="end">Iterator pointing past the last position vector of type Eigen::Vector3</param>
        /// <return>Second-order polynomial</return>
        template<typename forward_iterator_t>
        polynomial<typename std::iterator_traits<forward_iterator_t>::value_type::value_type> least_squares(forward_iterator_t begin, forward_iterator_t end);

        /// <summary>
        /// Solve least squares for Ax = b, returning x
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="rows">Number of rows of the matrix and entries of the vector</template>
        /// <template name="columns">Number of columns of the matrix and entries of the solution</template>
        /// <param name="matrix">Left-hand-side matrix</param>
        /// <param name="vector">Right-hand-side vector</param>
        /// <return>Solution, left-hand-side vector</return>
        template<typename floatp_t, int rows, int columns>
        Eigen::Matrix<floatp_t, columns, 1> least_squares(const Eigen::Matrix<floatp_t, rows, columns>& A, const Eigen::Matrix<floatp_t, rows, 1>& b);
    }
}

#include "tpf_least_squares.inl"
