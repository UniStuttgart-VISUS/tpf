#pragma once

#include "Eigen/Dense"

#include <functional>

namespace std
{
    /// <summary>
    /// Hash functor for coordinates
    /// </summary>
    template <>
    struct hash<Eigen::Matrix<long long, 1, 1>>
    {
        typedef Eigen::Matrix<long long, 1, 1> argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& coords) const;
    };

    /// <summary>
    /// Hash functor for coordinates
    /// </summary>
    template <>
    struct hash<Eigen::Matrix<long long, 2, 1>>
    {
        typedef Eigen::Matrix<long long, 2, 1> argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& coords) const;
    };

    /// <summary>
    /// Hash functor for coordinates
    /// </summary>
    template <>
    struct hash<Eigen::Matrix<long long, 3, 1>>
    {
        typedef Eigen::Matrix<long long, 3, 1> argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& coords) const;
    };
}

#include "tpf_hash.inl"