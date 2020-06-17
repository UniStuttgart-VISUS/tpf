#include "tpf_mpi_grid.h"

#include "tpf_mpi.h"
#include "tpf_mpi_exceptions.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_grid.h"
#include "../data/tpf_grid_information.h"

#include <algorithm>
#include <map>
#include <list>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        template <typename value_t, typename point_t, int dimensions, int components>
        inline data::grid<value_t, point_t, dimensions, components> mpi_grid::all_gather(const data::grid<value_t, point_t, dimensions, components>& source, MPI_Comm comm)
        {
#ifdef __tpf_use_mpi
            const auto& mpi = get_instance(comm);

            // Sanity checks
            if (!mpi.check_mpi_status() || source.get_size() == 0)
            {
                return source;
            }

            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(components > 0, "Number of components must be larger than zero");

            // Test for exceptions before continuing
            if (int errcode = mpi.communicate_error() != -1)
            {
                throw mpi_abort(errcode);
            }

            // Create target grid
            const auto global_extent = mpi_base::get_global_extent(source.get_extent(), comm);
            data::grid<value_t, point_t, dimensions, components> target(source.get_name(), global_extent);

            if (source.has_grid_information())
            {
                typename data::grid_information<point_t>::array_type cell_coordinates(dimensions);
                typename data::grid_information<point_t>::array_type node_coordinates(dimensions);
                typename data::grid_information<point_t>::array_type cell_sizes(dimensions);

                for (std::size_t d = 0; d < dimensions; ++d)
                {
                    const std::size_t num_elements = global_extent[d].second - global_extent[d].first + 1;

                    cell_coordinates[d].resize(num_elements);
                    node_coordinates[d].resize(num_elements + 1);
                    cell_sizes[d].resize(num_elements);
                }

                target.set_grid_information(std::move(cell_coordinates), std::move(node_coordinates), std::move(cell_sizes));
            }

            // Get local extent as subextent and extract subgrid from source
            auto subgrid = source.template extract_subgrid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type>(mpi_base::get_local_extent(source.get_extent(), comm).first);

            // In turn send and receive data
            for (int proc = 0; proc < mpi.get_num_processes(); ++proc)
            {
                if (proc == mpi.get_rank())
                {
                    // Broadcast
                    send_grid(subgrid, proc, true, source.has_grid_information(), MPI_COMM_WORLD);
                    target.fill_subgrid(subgrid, source.has_grid_information());
                }
                else
                {
                    // Receive broadcast
                    target.fill_subgrid(receive_grid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type, dimensions, components>
                        (proc, true, source.has_grid_information(), MPI_COMM_WORLD), source.has_grid_information());
                }
            }

            return target;
#else
            return source;
#endif
        }

        template <typename value_t, typename point_t, int dimensions, int components>
        inline void mpi_grid::synchronize_boundaries(data::grid<value_t, point_t, dimensions, components>& grid, MPI_Comm comm)
        {
#ifdef __tpf_use_mpi
            const auto& mpi = get_instance(comm);

            // Sanity checks
            if (!mpi.check_mpi_status())
            {
                return;
            }

            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(components > 0, "Number of components must be larger than zero");

            // Test for exceptions before continuing
            if (int errcode = mpi.communicate_error() != -1)
            {
                throw mpi_abort(errcode);
            }

            // Create topology if necessary
            const auto topology = create_topology(grid.get_extent(), comm);

            if (!topology.first)
            {
                return;
            }

            // Synchronize boundaries between two nodes in an independent set
            const int rank = mpi.get_rank();

            for (const auto& independent_set : topology.second.first)
            {
                for (const auto& ranks : independent_set)
                {
                    // Exchange boundaries
                    if (ranks.first == rank)
                    {
                        // Send, then receive
                        send_grid(grid.template extract_subgrid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type>(
                            topology.second.second.at(ranks.second)), ranks.second, false, false, comm);
                        grid.fill_subgrid(receive_grid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type, dimensions, components>(ranks.second, false, false, comm));
                    }
                    else if (ranks.second == rank)
                    {
                        // Receive, then send
                        grid.fill_subgrid(receive_grid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type, dimensions, components>(ranks.first, false, false, comm));
                        send_grid(grid.template extract_subgrid<typename mpi_t<value_t>::type, typename mpi_t<point_t>::type>(
                            topology.second.second.at(ranks.first)), ranks.first, false, false, comm);
                    }
                }

                mpi.barrier();
            }
#endif
        }

#ifdef __tpf_use_mpi
        inline std::pair<bool, std::pair<std::vector<std::vector<std::pair<int, int>>>, std::map<int, data::extent_t>>>
            mpi_grid::create_topology(const data::extent_t& extent, const MPI_Comm comm)
        {
            const auto& mpi = get_instance(comm);

            if (!mpi.check_mpi_status())
            {
                return std::make_pair(false, std::pair<std::vector<std::vector<std::pair<int, int>>>, std::map<int, data::extent_t>>());
            }

            // Clear previously set topology
            std::vector<std::vector<std::pair<int, int>>> topology;
            std::map<int, data::extent_t> border_extents;

            // Get all local extents
            const auto extents = mpi_base::get_extents(extent, comm);
            const auto local_extent = mpi_base::get_local_extent(extent, comm);

            if (local_extent.second == 0)
            {
                return std::make_pair(false, std::pair<std::vector<std::vector<std::pair<int, int>>>, std::map<int, data::extent_t>>());
            }

            // Find all neighboring processes, identified by their extent
            std::list<std::pair<int, int>> neighbors;

            const int own_rank = mpi.get_rank();

            for (int i = 0; i < extents.size() - 1; ++i)
            {
                for (int j = i + 1; j < extents.size(); ++j)
                {
                    // Do extents partially overlap?
                    bool overlap = true;

                    for (std::size_t dim = 0; dim < extent.size(); ++dim)
                    {
                        overlap &= extents[i][dim].first <= extents[j][dim].second && extents[i][dim].second >= extents[j][dim].first;
                    }

                    // Add them to the list
                    if (overlap)
                    {
                        neighbors.push_back(std::make_pair(i, j));

                        // If this thread is involved, calculate extent of shared boundary
                        const int remote_rank = (own_rank == i) ? j : i;

                        if (own_rank == i || own_rank == j)
                        {
                            // Calculate intersection of extents
                            data::extent_t intersection(extent.size());

                            for (std::size_t dim = 0; dim < extent.size(); ++dim)
                            {
                                intersection[dim].first = std::max(extent[dim].first, extents[remote_rank][dim].first);
                                intersection[dim].second = std::min(extent[dim].second, extents[remote_rank][dim].second);
                            }

                            // Get own boundary extent
                            data::extent_t boundary(extent.size());

                            for (std::size_t dim = 0; dim < extent.size(); ++dim)
                            {
                                boundary[dim].first = std::max(intersection[dim].first, local_extent.first[dim].first);
                                boundary[dim].second = std::min(intersection[dim].second, local_extent.first[dim].second);
                            }

                            // Add to border extents
                            border_extents.insert(std::make_pair(remote_rank, boundary));
                        }
                    }
                }
            }

            // Extract independent sets
            while (!neighbors.empty())
            {
                std::vector<std::pair<int, int>> independent_set;
                std::set<int> already_used;

                for (auto neighbor = neighbors.begin(); neighbor != neighbors.end();)
                {
                    if (already_used.find(neighbor->first) == already_used.find(neighbor->second))
                    {
                        independent_set.push_back(*neighbor);

                        already_used.insert(neighbor->first);
                        already_used.insert(neighbor->second);

                        neighbors.erase(neighbor++);
                    }
                    else
                    {
                        ++neighbor;
                    }
                }

                topology.push_back(independent_set);
            }

            return std::make_pair(true, std::make_pair(topology, border_extents));
        }

        template <typename value_t, typename point_t, int dimensions, int components>
        inline void mpi_grid::send_grid(const data::grid<value_t, point_t, dimensions, components>& source,
            const int receiver, const bool broadcast, const bool grid_info, const MPI_Comm comm)
        {
            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(components > 0, "Number of components must be larger than zero");
            static_assert(std::is_same<value_t, typename mpi_t<value_t>::type>::value, "Value type must be MPI compatible");
            static_assert(std::is_same<point_t, typename mpi_t<point_t>::type>::value, "Point type must be MPI compatible");

            const auto& mpi = get_instance(comm);

            // Send extent
            std::vector<typename mpi_t<std::size_t>::type> extent;

            for (const auto& ext : source.get_extent())
            {
                extent.push_back(static_cast<typename mpi_t<std::size_t>::type>(ext.first));
                extent.push_back(static_cast<typename mpi_t<std::size_t>::type>(ext.second));
            }

            if (broadcast)
            {
                mpi.broadcast(extent, mpi.get_rank());
            }
            else
            {
                mpi.send(extent, receiver);
            }

            // Send data
            if (broadcast)
            {
                mpi.broadcast(source.get_data(), mpi.get_rank());
            }
            else
            {
                mpi.send(source.get_data(), receiver);
            }

            // Send grid information
            if (grid_info)
            {
                for (std::size_t d = 0; d < dimensions; ++d)
                {
                    if (broadcast)
                    {
                        mpi.broadcast(source.get_cell_coordinates()[d], mpi.get_rank());
                        mpi.broadcast(source.get_node_coordinates()[d], mpi.get_rank());
                        mpi.broadcast(source.get_cell_sizes()[d], mpi.get_rank());
                    }
                    else
                    {
                        mpi.send(source.get_cell_coordinates()[d], receiver);
                        mpi.send(source.get_node_coordinates()[d], receiver);
                        mpi.send(source.get_cell_sizes()[d], receiver);
                    }
                }
            }
        }

        template <typename value_t, typename point_t, int dimensions, int components>
        inline data::grid<value_t, point_t, dimensions, components> mpi_grid::receive_grid(const int sender,
            const bool broadcast, const bool grid_info, const MPI_Comm comm)
        {
            static_assert(dimensions > 0, "Number of dimensions must be larger than zero");
            static_assert(components > 0, "Number of components must be larger than zero");
            static_assert(std::is_same<value_t, typename mpi_t<value_t>::type>::value, "Value type must be MPI compatible");
            static_assert(std::is_same<point_t, typename mpi_t<point_t>::type>::value, "Point type must be MPI compatible");

            const auto& mpi = get_instance(comm);

            // Receive extent
            std::vector<typename mpi_t<std::size_t>::type> extent_buffer;

            if (broadcast)
            {
                mpi.broadcast(extent_buffer, sender);
            }
            else
            {
                mpi.receive(extent_buffer, sender);
            }

            data::extent_t extent(dimensions);
            std::size_t index = 0;

            for (auto& ext : extent)
            {
                ext.first = static_cast<std::size_t>(extent_buffer[index++]);
                ext.second = static_cast<std::size_t>(extent_buffer[index++]);
            }

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

            data::grid<value_t, point_t, dimensions, components> received("Received", extent, std::move(data));

            // Receive grid information
            if (grid_info)
            {
                typename data::grid_information<point_t>::array_type cell_coordinates(dimensions);
                typename data::grid_information<point_t>::array_type node_coordinates(dimensions);
                typename data::grid_information<point_t>::array_type cell_sizes(dimensions);

                for (std::size_t d = 0; d < dimensions; ++d)
                {
                    if (broadcast)
                    {
                        mpi.broadcast(cell_coordinates[d], sender);
                        mpi.broadcast(node_coordinates[d], sender);
                        mpi.broadcast(cell_sizes[d], sender);
                    }
                    else
                    {
                        mpi.receive(cell_coordinates[d], sender);
                        mpi.receive(node_coordinates[d], sender);
                        mpi.receive(cell_sizes[d], sender);
                    }
                }

                received.set_grid_information(std::move(cell_coordinates), std::move(node_coordinates), std::move(cell_sizes));
            }

            return received;
        }
#endif
    }
}
