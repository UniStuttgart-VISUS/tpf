#include "tpf_hash.h"

#include "Eigen/Dense"

#include <limits>

namespace std
{
    inline typename hash<Eigen::Matrix<long long, 1, 1>>::result_type hash<Eigen::Matrix<long long, 1, 1>>::operator()(const argument_type& coords) const
    {
        return coords[0];
    }

    inline typename hash<Eigen::Matrix<long long, 2, 1>>::result_type hash<Eigen::Matrix<long long, 2, 1>>::operator()(const argument_type& coords) const
    {
        static const result_type stride = std::numeric_limits<result_type>::max() / 2;

        std::size_t hash = static_cast<std::size_t>(0);

        for (std::size_t d = 0; d < 2; ++d)
        {
            hash += d * stride + coords[d];
        }

        return hash;
    }

    inline typename hash<Eigen::Matrix<long long, 3, 1>>::result_type hash<Eigen::Matrix<long long, 3, 1>>::operator()(const argument_type& coords) const
    {
        static const result_type stride = std::numeric_limits<result_type>::max() / 3;

        std::size_t hash = static_cast<std::size_t>(0);

        for (std::size_t d = 0; d < 3; ++d)
        {
            hash += d * stride + coords[d];
        }

        return hash;
    }
}