#include "tpf_mpi_base.h"

#include "tpf_mpi.h"
#include "tpf_mpi_exceptions.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_data_information.h"

#include <cmath>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        inline std::pair<data::extent_t, std::size_t> mpi_base::get_local_extent(const data::extent_t& extent, const MPI_Comm comm)
        {
#ifdef __tpf_use_mpi
            const auto& mpi = get_instance(comm);

            if (!mpi.check_mpi_status())
            {
                return std::make_pair(extent, static_cast<std::size_t>(0));
            }

            // Test for exceptions before continuing
            if (int errcode = mpi.communicate_error() != -1)
            {
                throw mpi_abort(errcode);
            }

            // Get extents
            const auto global_extent = get_global_extent(extent, comm);
            auto extents = get_extents(extent, comm);

            // Increase number of ghost level count while decreasing extents until no overlap
            const auto num_processes = static_cast<std::size_t>(mpi.get_num_processes());
            const auto rank = static_cast<std::size_t>(mpi.get_rank());

            bool has_overlap = false;
            std::size_t num_ghost_levels = 0;

            do
            {
                // Find overlap
                has_overlap = false;

                for (std::size_t p = 0; p < num_processes && !has_overlap; ++p)
                {
                    if (rank != p)
                    {
                        bool has_overlap_local = true;

                        for (std::size_t d = 0; d < extent.size() && !has_overlap; ++d)
                        {
                            has_overlap_local &= extents[rank][d].first <= extents[p][d].second && extents[p][d].first <= extents[rank][d].second;
                        }

                        has_overlap |= has_overlap_local;
                    }
                }

                if (has_overlap)
                {
                    // Reduce extents
                    for (std::size_t p = 0; p < num_processes; ++p)
                    {
                        for (std::size_t d = 0; d < extent.size(); ++d)
                        {
                            if (extents[p][d].first > global_extent[d].first) ++extents[p][d].first;
                            if (extents[p][d].second < global_extent[d].second) --extents[p][d].second;
                        }
                    }

                    // Increase number of ghost levels
                    ++num_ghost_levels;
                }
            } while (has_overlap);

            return std::make_pair(extents[rank], num_ghost_levels);
#else
            return std::make_pair(extent, static_cast<std::size_t>(0));
#endif
        }

        inline std::vector<data::extent_t> mpi_base::split_space(const data::extent_t& extent, const std::size_t num_splits, const std::size_t num_ghost_levels)
        {
            // Recursive split on N processes:
            //  - Split room into "left" with floor(N / 2), and "right" with ceil(N / 2) processes
            //  - Assign width for "left" of extent * [N / floor(N / 2)], and for "right" extent * [N / ceil(N / 2)]
            //  - Recusively go on until number of split is one

            if (num_splits == 1)
            {
                std::vector<data::extent_t> output(1);
                output[0] = extent;

                return output;
            }

            const std::size_t x_num = extent[0].second - extent[0].first + 1;
            const std::size_t y_num = extent[1].second - extent[1].first + 1;
            const std::size_t z_num = extent[2].second - extent[2].first + 1;

            const std::size_t largest_direction = (x_num > y_num && x_num > z_num) ? 0 : ((y_num > z_num) ? 1 : 2);

            const std::size_t num_left = static_cast<std::size_t>(std::floor(static_cast<float>(num_splits) / 2.0f));
            const std::size_t num_right = static_cast<std::size_t>(std::ceil(static_cast<float>(num_splits) / 2.0f));

            data::extent_t sub_left = extent;
            data::extent_t sub_right = extent;

            sub_left[largest_direction].second = sub_left[largest_direction].first +
                static_cast<std::size_t>(static_cast<float>((sub_left[largest_direction].second
                    - sub_left[largest_direction].first)) * static_cast<float>(num_left) / static_cast<float>(num_splits));
            sub_right[largest_direction].first = sub_left[largest_direction].second + 1;

            sub_left[largest_direction].second += num_ghost_levels;
            sub_right[largest_direction].first -= num_ghost_levels;

            auto left = split_space(sub_left, num_left, num_ghost_levels);
            auto right = split_space(sub_right, num_right, num_ghost_levels);

            left.insert(left.end(), right.begin(), right.end());

            return left;
        }

#ifdef __tpf_use_mpi
        inline data::extent_t mpi_base::get_global_extent(const data::extent_t& extent, const MPI_Comm comm)
        {
            const auto& mpi = get_instance(comm);

            // If MPI is not active, return extent
            data::extent_t global_extent = extent;

            if (!mpi.check_mpi_status())
            {
                return global_extent;
            }

            // Reduce all to global extent
            for (std::size_t d = 0; d < extent.size(); ++d)
            {
                mpi.allreduce_inplace(global_extent[d].first, MPI_MIN);
                mpi.allreduce_inplace(global_extent[d].second, MPI_MAX);
            }

            return global_extent;
        }

        inline std::vector<data::extent_t> mpi_base::get_extents(const data::extent_t& extent, const MPI_Comm comm)
        {
            const auto& mpi = get_instance(comm);

            // If MPI is not active, return extent
            std::vector<data::extent_t> local_extents;
            local_extents.push_back(extent);

            if (!mpi.check_mpi_status())
            {
                return local_extents;
            }

            // Copy extent into contiguous memory
            std::vector<typename mpi_t<std::size_t>::type> extent_copy;

            for (auto& ext : extent)
            {
                extent_copy.push_back(static_cast<typename mpi_t<std::size_t>::type>(ext.first));
                extent_copy.push_back(static_cast<typename mpi_t<std::size_t>::type>(ext.second));
            }

            // Gather all
            std::vector<typename mpi_t<std::size_t>::type> local_extents_copy;

            mpi.allgather(extent_copy, local_extents_copy);

            // Unpack from receive buffer
            local_extents.resize(static_cast<std::size_t>(mpi.get_num_processes()));
            std::size_t i = 0;

            for (auto& local_extent : local_extents)
            {
                local_extent.resize(extent.size());

                for (auto& ext : local_extent)
                {
                    ext.first = static_cast<std::size_t>(local_extents_copy[i++]);
                    ext.second = static_cast<std::size_t>(local_extents_copy[i++]);
                }
            }

            return local_extents;
        }

        template <typename value_t>
        inline void mpi_base::send_array(const std::vector<value_t>& source, const int receiver, const bool broadcast, const MPI_Comm comm)
        {
            const auto& mpi = get_instance(comm);

            // Send data
            if (broadcast)
            {
                mpi.broadcast(source, mpi.get_rank());
            }
            else
            {
                mpi.send(source, receiver);
            }
        }

        template <typename value_t>
        inline std::vector<value_t> mpi_base::receive_array(const int sender, const bool broadcast, const MPI_Comm comm)
        {
            const auto& mpi = get_instance(comm);

            // Receive data
            std::vector<value_t> data;

            if (broadcast)
            {
                mpi.broadcast(data, sender);
            }
            else
            {
                mpi.receive(data, sender);
            }

            // Set data
            return data;
        }
#endif
    }
}
