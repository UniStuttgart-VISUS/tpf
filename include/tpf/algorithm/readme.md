# Algorithm

Here, general algorithms are provided that are independent from the TPF data structures.

- [Jenkins one at a time hash](#jenkins-one-at-a-time-hash)
- [Least squares](#least-squares)
- [Loop](#loop)

## Jenkins one at a time hash

This function computes a combined hash from the input values using the [Jenkins hash function](https://en.wikipedia.org/wiki/Jenkins_hash_function).

```c++
template <typename first_t, typename second_t = void*>
std::uint32_t
joaat_hash(const first_t& first, const second_t second = nullptr);

template <typename first_t, typename second_t, typename... other_t>
std::uint32_t
joaat_hash(const first_t& first, const second_t& second, const other_t&... values);
```

The first function computes the hash for one or two input values, the second allows a user-defined number of parameters in a variadic template.

## Least squares

This approximates the least squares for a second-order polynomial.

### polynomial

```c++
template <typename floatp_t>
struct polynomial
{
    floatp_t a_xx;
    floatp_t a_yy;
    floatp_t a_xy;
    floatp_t a_x;
    floatp_t a_y;
    floatp_t a_1;

    polynomial<floatp_t> operator-(const polynomial<floatp_t>& poly) const;
};
```

This structure stores the coefficients of the second-order polynomial of the form
``a_xx`` * x² + ``a_yy`` * y² + ``a_xy`` * x * y + ``a_x`` * x + ``a_y`` * y + ``a_1``.

Additionally, it overloads the ``-`` operator for subtraction of polynomials.

### least_squares

As input, an iterator range for a collection of ``Eigen::Vector3`` is expected, for which a polynomial is computed as least squares approximation.

```c++
template<typename forward_iterator_t>
polynomial<typename std::iterator_traits<forward_iterator_t>::value_type::value_type>
least_squares(forward_iterator_t begin, forward_iterator_t end);
```

In the second function, the least squares approximation is performed for the equation ``Ax = b``, solving for ``x``.

```c++
template<typename floatp_t, int rows, int columns>
Eigen::Matrix<floatp_t, columns, 1>
least_squares(const Eigen::Matrix<floatp_t, rows, columns>& A, const Eigen::Matrix<floatp_t, rows, 1>& b);
```

## Loop

This implements a nested loop over all elements defined by the extents given by ``start`` and ``end``. The number of nested loops corresponds to the dimension of these extents. For each element, the provided ``function`` is called with its calculated index as parameter. The result of the function is written using the output iterator ``out_begin``, if provided. The parallel version uses OpenMP for parallelization of the outer loop, with a critical section for the output, ensuring thread-safety.

```c++
template <typename integral_t,
          std::size_t dimensions,
          typename function_t,
          typename output_iterator_t>
void
nested_loop(const Eigen::Matrix<integral_t, dimensions, 1>& start,
            const Eigen::Matrix<integral_t, dimensions, 1>& end,
            function_t function,
            output_iterator_t& out_begin = output_iterator_t());

template <typename integral_t,
          std::size_t dimensions,
          typename function_t,
          typename output_iterator_t>
void
parallel_nested_loop(
              const Eigen::Matrix<integral_t, dimensions, 1>& start,
              const Eigen::Matrix<integral_t, dimensions, 1>& end,
              function_t function,
              output_iterator_t& out_begin = output_iterator_t());
```

