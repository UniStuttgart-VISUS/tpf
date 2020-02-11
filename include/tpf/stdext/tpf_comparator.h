#pragma once

#include "../data/tpf_data_information.h"

#include "Eigen/Dense"

#include <functional>

namespace std
{
    /// <summary>
    /// Less functor for extents
    /// </summary>
    template <>
    struct less<tpf::data::extent_t>
    {
        typedef bool result_type;
        typedef tpf::data::extent_t first_argument_type;
        typedef first_argument_type second_argument_type;

        bool operator()(const first_argument_type& lhs, const second_argument_type& rhs) const;
    };

    /// <summary>
    /// Less functor for coordinates
    /// </summary>
    template <>
    struct less<Eigen::Matrix<long long, 1, 1>>
    {
        typedef bool result_type;
        typedef Eigen::Matrix<long long, 1, 1> first_argument_type;
        typedef first_argument_type second_argument_type;

        bool operator()(const first_argument_type& lhs, const second_argument_type& rhs) const;
    };

    /// <summary>
    /// Less functor for coordinates
    /// </summary>
    template <>
    struct less<Eigen::Matrix<long long, 2, 1>>
    {
        typedef bool result_type;
        typedef Eigen::Matrix<long long, 2, 1> first_argument_type;
        typedef first_argument_type second_argument_type;

        bool operator()(const first_argument_type& lhs, const second_argument_type& rhs) const;
    };

    /// <summary>
    /// Less functor for coordinates
    /// </summary>
    template <>
    struct less<Eigen::Matrix<long long, 3, 1>>
    {
        typedef bool result_type;
        typedef Eigen::Matrix<long long, 3, 1> first_argument_type;
        typedef first_argument_type second_argument_type;

        bool operator()(const first_argument_type& lhs, const second_argument_type& rhs) const;
    };

    /// <summary>
    /// Less functor for points
    /// Comparision of points is not physical meaningful, but useful to make sets and maps of points.
    /// Therefore the order of the points doesn't really matter. One easy compare function is first order by x-dimension, then y-dim, then z-dim.
    /// </summary>
    template <typename floatp_t>
    struct less<Eigen::Matrix<floatp_t, 3, 1>>
    {
        typedef bool result_type;
        typedef Eigen::Matrix<floatp_t, 3, 1> first_argument_type;
        typedef first_argument_type second_argument_type;

        bool operator()(const first_argument_type& lhs, const second_argument_type& rhs) const;
    };
}

#include "tpf_comparator.inl"