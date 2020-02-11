#pragma once

#include "tpf_mpi.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_data_information.h"

#include <utility>
#include <map>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        /// <summary>
        /// MPI class for basic high-level functions
        /// </summary>
        class mpi_base
        {
        public:
            /// <summary>
            /// Calculate the local extent, excluding synchronizable boundaries
            /// </summary>
            /// <param name="extent">Extent</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Local extent and number of ghost levels</returns>
            static std::pair<data::extent_t, std::size_t> get_local_extent(const data::extent_t& extent, MPI_Comm comm = MPI_COMM_WORLD);

            /// <summary>
            /// Split space into subspaces
            /// </summary>
            /// <param name="extent">Whole extent</param>
            /// <param name="num_splits">Number of splits</param>
            /// <param name="num_ghost_levels">Number of ghost levels</param>
            /// <returns>Extents for each process including ghost cells</returns>
            static std::vector<data::extent_t> split_space(const data::extent_t& extent, std::size_t num_splits, std::size_t num_ghost_levels);

#ifdef __tpf_use_mpi
        protected:
            /// <summary>
            /// Calculate the global extent
            /// </summary>
            /// <param name="extent">Local extent</param>
            /// <param name="comm">Communicator</param>
            /// <return>Global extent</return>
            static data::extent_t get_global_extent(const data::extent_t& extent, MPI_Comm comm);

            /// <summary>
            /// Gather extents
            /// </summary>
            /// <param name="extent">Extent</param>
            /// <param name="comm">Communicator</param>
            /// <return>Extents</return>
            static std::vector<data::extent_t> get_extents(const data::extent_t& extent, MPI_Comm comm);

            /// <summary>
            /// Send an array
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <param name="source">Source array</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            template <typename value_t>
            static void send_array(const std::vector<value_t>& source, int receiver, bool broadcast, MPI_Comm comm);

            /// <summary>
            /// Receive (a part of) an array
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <param name="sender">Rank of the sender</param>
            /// <param name="broadcast">True: receive from broadcast, false: receive from send</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Received array</returns>
            template <typename value_t>
            static std::vector<value_t> receive_array(int sender, bool broadcast, MPI_Comm comm);
#endif
        };
    }
}

#include "tpf_mpi_base.inl"