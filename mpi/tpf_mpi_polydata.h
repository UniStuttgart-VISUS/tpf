#pragma once

#include "tpf_mpi_base.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_polydata.h"
#include "../data/tpf_data_information.h"

#include "../geometry/tpf_geometric_object.h"

#include <memory>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        /// <summary>
        /// MPI class for polydata-related high-level functions
        /// </summary>
        class mpi_polydata : public mpi_base
        {
        public:
            /// <summary>
            /// MPI allgather for poly data objects
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="KernelPT">CGAL kernel type</template>
            /// <template name="data_information_t">Data information type</template>
            /// <param name="source">Source poly data</param>
            /// <param name="comm">Communicator</param>
            /// <param name="data_infos">Data information for arrays to gather</param>
            /// <return>Gathered poly data</return>
            template <typename point_t, typename kernel_t = geometry::default_kernel_t, typename... data_information_t>
            static data::polydata<point_t> all_gather(const data::polydata<point_t>& source, MPI_Comm comm = MPI_COMM_WORLD, const data_information_t&... data_infos);

#ifdef __tpf_use_mpi
        private:
            /// <summary>
            /// Send poly data
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <param name="source">Source poly data</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            template <typename point_t>
            static void send_polydata(const data::polydata<point_t>& source, int receiver, bool broadcast, MPI_Comm comm);

            /// <summary>
            /// Send poly data
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="data_information_t">Data information objects</template>
            /// <param name="source">Source poly data</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="data_infos">Data information objects</param>
            template <typename point_t, typename... data_information_t>
            static void send_polydata(const data::polydata<point_t>& source, int receiver,
                bool broadcast, MPI_Comm comm, const data_information_t&... data_infos);

            /// <summary>
            /// Send a geometric object
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <param name="source">Source object</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            template <typename point_t>
            static void send_geometry(const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& source,
                int receiver, bool broadcast, MPI_Comm comm);

            /// <summary>
            /// Send selected poly data arrays
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="value_t">Value type of the attached array</template>
            /// <template name="rows">Number of rows of the attached array</template>
            /// <template name="columns">Number of columns of the attached array</template>
            /// <template name="data_information_t">Data information objects</template>
            /// <param name="source">Source poly data</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="first_data_info">First data information object</param>
            /// <param name="more_data_infos">More data information objects</param>
            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_information_t>
            static void send_polydata_arrays(const data::polydata<point_t>& source, int receiver, bool broadcast,
                MPI_Comm comm, const data::data_information<value_t, rows, columns>& first_data_info, const data_information_t&... more_data_infos);

            /// <summary>
            /// Send selected poly data arrays
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="value_t">Value type of the attached array</template>
            /// <template name="rows">Number of rows of the attached array</template>
            /// <template name="columns">Number of columns of the attached array</template>
            /// <param name="source">Source poly data</param>
            /// <param name="receiver">Rank of the receiver (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="first_data_info">Data information object</param>
            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns>
            static void send_polydata_arrays(const data::polydata<point_t>& source, int receiver, bool broadcast,
                MPI_Comm comm, const data::data_information<value_t, rows, columns>& data_info);

            /// <summary>
            /// Receive poly data
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="KernelPT">CGAL kernel type</template>
            /// <param name="sender">Rank of the sender</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Received poly data</returns>
            template <typename point_t, typename kernel_t>
            static data::polydata<point_t> receive_polydata(int sender, bool broadcast, MPI_Comm comm);

            /// <summary>
            /// Receive poly data
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="KernelPT">CGAL kernel type</template>
            /// <template name="data_information_t">Data information objects</template>
            /// <param name="sender">Rank of the sender</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="data_infos">Data information objects</param>
            /// <returns>Received poly data</returns>
            template <typename point_t, typename kernel_t, typename... data_information_t>
            static data::polydata<point_t> receive_polydata(int sender, bool broadcast, MPI_Comm comm, const data_information_t&... data_infos);

            /// <summary>
            /// Receive a geometric object
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="KernelPT">CGAL kernel type</template>
            /// <param name="sender">Rank of the sender</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <returns>Received geometric object</returns>
            template <typename point_t, typename kernel_t>
            static std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> receive_geometry(int sender, bool broadcast, MPI_Comm comm);

            /// <summary>
            /// Receive selected poly data arrays
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="value_t">Value type of the attached array</template>
            /// <template name="rows">Number of rows of the attached array</template>
            /// <template name="columns">Number of columns of the attached array</template>
            /// <template name="data_information_t">Data information objects</template>
            /// <param name="target">Target poly data</param>
            /// <param name="sender">Rank of the sender (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="first_data_info">First data information object</param>
            /// <param name="more_data_infos">More data information objects</param>
            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_information_t>
            static void receive_polydata_arrays(data::polydata<point_t>& target, int sender, bool broadcast,
                MPI_Comm comm, const data::data_information<value_t, rows, columns>& first_data_info, const data_information_t&... more_data_infos);

            /// <summary>
            /// Receive selected poly data arrays
            /// </summary>
            /// <template name="point_t">Point type</template>
            /// <template name="value_t">Value type of the attached array</template>
            /// <template name="rows">Number of rows of the attached array</template>
            /// <template name="columns">Number of columns of the attached array</template>
            /// <param name="target">Target poly data</param>
            /// <param name="sender">Rank of the sender (ignored when broadcasting)</param>
            /// <param name="broadcast">True: send as broadcast, false: send to given receiver</param>
            /// <param name="comm">Communicator</param>
            /// <param name="first_data_info">Data information object</param>
            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns>
            static void receive_polydata_arrays(data::polydata<point_t>& target, int sender, bool broadcast,
                MPI_Comm comm, const data::data_information<value_t, rows, columns>& data_info);
#endif
        };
    }
}

#include "tpf_mpi_polydata.inl"