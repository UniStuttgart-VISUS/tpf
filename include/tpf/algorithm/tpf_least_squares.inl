#include "tpf_least_squares.h"

#include "../log/tpf_log.h"

#include "Eigen/Dense"

#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace tpf
{
    namespace algorithm
    {
        template <typename floatp_t>
        inline polynomial<floatp_t> polynomial<floatp_t>::operator-(const polynomial<floatp_t>& poly) const
        {
            polynomial<floatp_t> new_poly;

            new_poly.a_xx = a_xx - poly.a_xx;
            new_poly.a_yy = a_yy - poly.a_yy;
            new_poly.a_xy = a_xy - poly.a_xy;
            new_poly.a_x = a_x - poly.a_x;
            new_poly.a_y = a_y - poly.a_y;
            new_poly.a_1 = a_1 - poly.a_1;

            return new_poly;
        }

        template<typename forward_iterator_t>
        inline polynomial<typename std::iterator_traits<forward_iterator_t>::value_type::value_type> least_squares(forward_iterator_t current, const forward_iterator_t end)
        {
            using floatp_t = typename std::iterator_traits<forward_iterator_t>::value_type::value_type;

            static_assert(std::is_floating_point<floatp_t>::value,
                "Type must be floating point");
            static_assert(std::is_same<typename std::iterator_traits<forward_iterator_t>::value_type, Eigen::Matrix<floatp_t, 3, 1>>::value,
                "The value type of the iterator does not match the requirements");
            static_assert(std::is_base_of<std::forward_iterator_tag, typename std::iterator_traits<forward_iterator_t>::iterator_category>::value,
                "Iterator must be a forward input iterator");
            
            // Check number of positions
            const auto num_positions = std::distance(current, end);

            if (num_positions < 6)
            {
                throw std::runtime_error(__tpf_error_message("Least squares needs at least 6 positions to compute a 2nd-order polynomial"));
            }

            // Create matrix A and right-hand side b
            Eigen::Matrix<floatp_t, Eigen::Dynamic, 6> A;
            A.resize(num_positions, Eigen::NoChange);

            Eigen::Matrix<floatp_t, Eigen::Dynamic, 1> b;
            b.resize(num_positions, Eigen::NoChange);

            for (std::size_t j = 0; current != end; ++current, ++j)
            {
                A(j, 0) = 1;
                A(j, 1) = (*current)[0];
                A(j, 2) = (*current)[1];
                A(j, 3) = (*current)[0] * (*current)[1];
                A(j, 4) = (*current)[0] * (*current)[0];
                A(j, 5) = (*current)[1] * (*current)[1];

                b(j) = (*current)[2];
            }

            // Solve least squares for Ax = b
            const auto result = least_squares(A, b);

            polynomial<floatp_t> poly;
            poly.a_1 = result[0];
            poly.a_x = result[1];
            poly.a_y = result[2];
            poly.a_xy = result[3];
            poly.a_xx = result[4];
            poly.a_yy = result[5];

            return poly;
        }

        template<typename floatp_t, int rows, int columns>
        inline Eigen::Matrix<floatp_t, columns, 1> least_squares(const Eigen::Matrix<floatp_t, rows, columns>& A, const Eigen::Matrix<floatp_t, rows, 1>& b)
        {
            // Solve least squares for Ax = b
            return Eigen::JacobiSVD<Eigen::Matrix<floatp_t, rows, columns>>(A, Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
        }
    }
}
