#pragma once

#include "tpf_mpi_base.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_data_information.h"
#include "../data/tpf_grid.h"

#include <map>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        /// <summary>
        /// MPI class for grid-related high-level functions
        /// </summary>
        class mpi_grid : public mpi_base
        {
        public:
            /// <summary>
            /// MPI allgather for arrays with potential ghost cells
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            /// <template name="dimensions">Number of spatial dimensions</template>
            /// <template name="components">Number of components</template>
            /// <param name="source">Source grid</param>
            /// <param name="comm">Communicator</param>
            /// <return>Gathered grid</return>
            template <typename value_t, typename point_t, int dimensions, int components>
            static data::grid<value_t, point_t, dimensions, components> all_gather(const data::grid<value_t, point_t, dimensions, components>& source, MPI_Comm comm = MPI_COMM_WORLD);

            /// <summary>
            /// Synchronize boundaries (ghost cells)
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            /// <template name="dimensions">Number of spatial dimensions</template>
            /// <template name="components">Number of components</template>
            /// <param name="source">Input/output data grid</param>
            /// <param name="comm">Communicator</param>
            template <typename value_t, typename point_t, int dimensions, int components>
            static void synchronize_boundaries(data::grid<value_t, point_t, dimensions, components>& source, MPI_Comm comm = MPI_COMM_WORLD);

#ifdef __tpf_use_mpi
        private:
            /// <summary>
            /// Create topology for communciation between neighbors
            /// </summary>
            /// <param name="extent">Extent</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Topology [success; [topology; border extents]]</returns>
            static std::pair<bool, std::pair<std::vector<std::vector<std::pair<int, int>>>, std::map<int, data::extent_t>>> create_topology(const data::extent_t& extent, MPI_Comm comm);

            /// <summary>
            /// Send a grid
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            /// <template name="dimensions">Number of spatial dimensions</template>
            /// <template name="components">Number of components</template>
            /// <param name="source">Source grid</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="grid_info">True: send also grid information, false: send grid data only</param>
            /// <param name="comm">Communicator</param>
            template <typename value_t, typename point_t, int dimensions, int components>
            static void send_grid(const data::grid<value_t, point_t, dimensions, components>& source, int receiver, bool broadcast, bool grid_info, MPI_Comm comm);

            /// <summary>
            /// Receive (part of) a grid
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="point_t">Point type</template>
            /// <template name="dimensions">Number of spatial dimensions</template>
            /// <template name="components">Number of components</template>
            /// <param name="sender">Rank of the sender</param>
            /// <param name="broadcast">True: receive from broadcast, false: receive from send</param>
            /// <param name="grid_info">True: receive also grid information, false: receive grid data only</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Received grid</returns>
            template <typename value_t, typename point_t, int dimensions, int components>
            static data::grid<value_t, point_t, dimensions, components> receive_grid(int sender, bool broadcast, bool grid_info, MPI_Comm comm);
#endif
        };
    }
}

#include "tpf_mpi_grid.inl"
