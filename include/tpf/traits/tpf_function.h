#pragma once

#include <type_traits>

/// Macro for creating a class to check if a certain function exists
#define __tpf_has_function(function) \
template <typename __T> \
class has_function_ ## function \
{ \
    template <typename C> static char test(decltype(&C::function)); \
    template <typename C> static long test(...); \
    \
public: \
    enum { value = sizeof(test<__T>(0)) == sizeof(char) }; \
};

/// Macro for creating a class to check if a certain operator exists
#define __tpf_has_operator(op, sign, name) \
struct __tpf_no_ ## name {}; \
template<typename lhs_t, typename rhs_t> __tpf_no_ ## name op (const lhs_t&, const rhs_t&); \
template<typename lhs_t, typename rhs_t = lhs_t> \
struct has_operator_ ## name \
{ \
    enum { value = !std::is_same<decltype(std::declval<lhs_t>() sign std::declval<rhs_t>()), __tpf_no_ ## name >::value }; \
};
