#include "tpf_comparator.h"

#include "../data/tpf_data_information.h"

#include "Eigen/Dense"

namespace std
{
    inline bool less<tpf::data::extent_t>::operator()(const first_argument_type& lhs, const second_argument_type& rhs) const
    {
        if (lhs.size() != rhs.size())
        {
            return false;
        }

        for (std::size_t d = 0; d < lhs.size(); ++d)
        {
            if (lhs[lhs.size() - d - 1].first < rhs[lhs.size() - d - 1].first)
            {
                return true;
            }
            else if (lhs[lhs.size() - d - 1].first > rhs[lhs.size() - d - 1].first)
            {
                return false;
            }
        }

        return false;
    }

    inline bool less<Eigen::Matrix<long long, 1, 1>>::operator()(const first_argument_type& lhs, const second_argument_type& rhs) const
    {
        return lhs[0] < rhs[0];
    }

    inline bool less<Eigen::Matrix<long long, 2, 1>>::operator()(const first_argument_type& lhs, const second_argument_type& rhs) const
    {
        return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
    }

    inline bool less<Eigen::Matrix<long long, 3, 1>>::operator()(const first_argument_type& lhs, const second_argument_type& rhs) const
    {
        return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]) || (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] < rhs[2]);
    }

    template <typename floatp_t>
    inline bool less<Eigen::Matrix<floatp_t, 3, 1>>::operator()(const first_argument_type& lhs, const second_argument_type& rhs) const
    {
        if (lhs[0] == rhs[0]) {
            if (lhs[1] == rhs[1]) {
                return lhs[2] < rhs[2];
            }
            return lhs[1] < rhs[1];
        }
        return lhs[0] < rhs[0];
    };
}