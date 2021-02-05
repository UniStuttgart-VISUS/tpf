# Utility

This provides utility functions:

- [Copy](#copy)
- [Downcast](#downcast)
- [Min-Max](#min-max)
- [Optional](#optional)
- [Value Cast](#value-cast)

## Copy

The ``copy `` function copies the values from one vector to the other, invoking their ``clone`` function for a deep copy.

```c++
template <typename value_t>
void
copy(const std::vector<std::shared_ptr<value_t>>& from, std::vector<std::shared_ptr<value_t>>& to);
```

The ``copy_data`` function copies data only if the data types are different, by casting each element. If the data types are the same, the data is returned by reference without modification.

```c++
template <typename new_t, typename original_t>
std::vector<new_t>&
copy_data(std::vector<original_t>& input,
          typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
const std::vector<new_t>&
copy_data(const std::vector<original_t>& input,
          typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
std::vector<new_t>
copy_data(const std::vector<original_t>& input,
          typename std::enable_if<!std::is_same<original_t, new_t>::value>::type* = nullptr);
```

``copy_grid`` behaves like ``copy_data``, but works on a nested vector of data, i.e., a 2D data field.

```c++
template <typename new_t, typename original_t>
std::vector<std::vector<new_t>>&
copy_grid(std::vector<std::vector<original_t>>& input,
          typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
const std::vector<std::vector<new_t>>&
copy_grid(const std::vector<std::vector<original_t>>& input,
          typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
std::vector<std::vector<new_t>>
copy_grid(const std::vector<std::vector<original_t>>& input,
          typename std::enable_if<!std::is_same<original_t, new_t>::value>::type* = nullptr);
```

## Downcast

This cast checks the data types, and performs sanity checks if the cast is narrowing. In this case, an exception is thrown if any of the following condition meets:

- the original value is negative, the new type is signed,
- the original value is greater than the maximum value of the new type, or
- the original value is smaller than the lowest value of the new type.

If the cast is to the same type or a widening cast, no checks are performed, but the original value returned or cast into the new type, respectively.

```c++
template <typename new_t, typename original_t>
new_t
safe_downcast(original_t original,
              typename std::enable_if<narrowing<new_t, original_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
new_t
safe_downcast(const original_t& original,
              typename std::enable_if<std::is_same<new_t, original_t>::value>::type* = nullptr);

template <typename new_t, typename original_t>
new_t
safe_downcast(const original_t& original,
              typename std::enable_if<widening<new_t, original_t>::value>::type* = nullptr);
```

## Min-Max

Variadic implementation of minimum and maximum functions for a varying number of input parameters.

```c++
template <typename value_t>
value_t
min(value_t first, value_t second);

template <typename value_t, typename... other_t>
value_t
min(value_t first, value_t second, other_t... others);

template <typename value_t>
value_t
max(value_t first, value_t second);

template <typename value_t, typename... other_t>
value_t
max(value_t first, value_t second, other_t... others);
```

Find the minimum or maximum in a range of values.
*Note that these functions are deprecated, as the C++ standard provides [``std::min_element``](https://en.cppreference.com/w/cpp/algorithm/min_element) and [``std::max_element``](https://en.cppreference.com/w/cpp/algorithm/max_element) with the same functionality.*

```c++
template <typename iterator_t, typename comparator = std::less<typename iterator_t::value_type>>
auto
min_in_range(iterator_t begin, iterator_t end, comparator comp = std::less<typename iterator_t::value_type>());

template <typename iterator_t, typename comparator = std::less<typename iterator_t::value_type>>
auto
max_in_range(iterator_t begin, iterator_t end, comparator comp = std::greater<typename iterator_t::value_type>());
```

## Optional

This class provides an implementation based on the idea of [``std::optional``](https://en.cppreference.com/w/cpp/utility/optional), which was introduced with C++ 17. This code can be used with older C++ versions and hence provides compatibility.

## Value Cast

This cast function uses ``static_cast`` on a value, or on each value of a vector. In case the types are identical, a constant reference or copy of the element is returned.

```c++
template <typename new_t, typename original_t>
constexpr new_t
value_cast(original_t original);

template <typename new_value_t, typename original_value_t, int rows, int columns>
constexpr Eigen::Matrix<new_value_t, rows, columns>
value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
           typename std::enable_if<!std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);

template <typename new_value_t, typename original_value_t, int rows, int columns>
constexpr const Eigen::Matrix<new_value_t, rows, columns>&
value_cast(const Eigen::Matrix<original_value_t, rows, columns>& original,
           typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);

template <typename new_value_t, typename original_value_t, int rows, int columns>
constexpr Eigen::Matrix<new_value_t, rows, columns>
value_cast(Eigen::Matrix<original_value_t, rows, columns> original,
           typename std::enable_if<std::is_same<original_value_t, new_value_t>::value>::type* = nullptr);
```

