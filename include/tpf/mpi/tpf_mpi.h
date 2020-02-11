#pragma once

#ifdef __tpf_use_mpi
#include "mpi.h"
#endif

#include "tpf_mpi_traits.h"

#include <map>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        /// <summary>
        /// MPI wrapper
        /// </summary>
        class mpi
        {
        public:
            /// <summary>
            /// Singleton instance
            /// </summary>
            /// <param name="comm">Communicator</param>
            /// <return>mpi instance</return>
            static const mpi& get_instance(MPI_Comm comm = MPI_COMM_WORLD);

            /// Deleted copy, move constructor and assignment, move operator
            mpi() = delete;
            mpi(const mpi&) = delete;
            mpi(mpi&&) = delete;
            mpi& operator=(const mpi&) = delete;
            mpi& operator=(mpi&&) = delete;

            /// <summary>
            /// Destructor
            /// </summary>
            ~mpi();

            /// <summary>
            /// Check for enabled mpi
            /// </summary>
            /// <return>True if mpi is enabled, false otherwise</return>
            bool check_mpi_status() const;

            /// <summary>
            /// Return the rank
            /// </summary>
            /// <returns>Rank</returns>
            int get_rank() const;

            /// <summary>
            /// Return the number of processes
            /// </summary>
            /// <returns>Number of processes</returns>
            int get_num_processes() const;

            /// <summary>
            /// Return the communicator
            /// </summary>
            /// <returns>Communicator</returns>
            MPI_Comm get_comm() const;

            /// <summary>
            /// Communicate a possible error
            /// </summary>
            /// <param name="error">Has there been an error?</param>
            /// <returns>Rank of a failed process; or -1 for no error</returns>
            int communicate_error(bool error = false) const;

            /// <summary>
            /// Broadcast single value
            /// </summary>
            /// <template name="value_t">Value type of the in- or output</template>
            /// <param name="value">Input or output value</param>
            /// <param name="root">Sending rank</param>
            template <typename value_t>
            void broadcast(value_t& value, int root) const;

            /// <summary>
            /// Broadcast single value (only send)
            /// </summary>
            /// <template name="value_t">Value type of the in- or output</template>
            /// <param name="value">Input or output value</param>
            /// <param name="root">Sending rank</param>
            template <typename value_t>
            void broadcast(const value_t& value, int root) const;

            /// <summary>
            /// Broadcast an array of values
            /// </summary>
            /// <template name="value_t">Value type of the in- or output</template>
            /// <param name="values">Input or output values</param>
            /// <param name="root">Sending rank</param>
            template <typename value_t>
            void broadcast(std::vector<value_t>& values, int root) const;

            /// <summary>
            /// Broadcast an array of values (only send)
            /// </summary>
            /// <template name="value_t">Value type of the in- or output</template>
            /// <param name="values">Input or output values</param>
            /// <param name="root">Sending rank</param>
            template <typename value_t>
            void broadcast(const std::vector<value_t>& values, int root) const;

            /// <summary>
            /// Send single value
            /// </summary>
            /// <template name="value_t">Value type of the input</template>
            /// <param name="value">Input value</param>
            /// <param name="destination">Receiving rank</param>
            template <typename value_t>
            void send(const value_t& value, int destination) const;

            /// <summary>
            /// Send an array of values
            /// </summary>
            /// <template name="value_t">Value type of the input</template>
            /// <param name="values">Input values</param>
            /// <param name="destination">Receiving rank</param>
            template <typename value_t>
            void send(const std::vector<value_t>& values, int destination) const;

            /// <summary>
            /// Receive single value
            /// </summary>
            /// <template name="value_t">Value type of the output</template>
            /// <param name="value">Output value</param>
            /// <param name="source">Sending rank</param>
            template <typename value_t>
            void receive(value_t& value, int source) const;

            /// <summary>
            /// Receive an array of values
            /// </summary>
            /// <template name="value_t">Value type of the output</template>
            /// <param name="values">Output values</param>
            /// <param name="source">Sending rank</param>
            template <typename value_t>
            void receive(std::vector<value_t>& values, int source) const;

            /// <summary>
            /// Allgather a value per process
            /// </summary>
            /// <template name="value_t">Value type of the in- and output</template>
            /// <param name="value">Input value</param>
            /// <param name="output">Gathered output values</param>
            template <typename value_t>
            void allgather(const value_t& value, std::vector<value_t>& output) const;

            /// <summary>
            /// Allgather multiple values per process
            /// </summary>
            /// <template name="value_t">Value type of the in- and output</template>
            /// <param name="values">Input values</param>
            /// <param name="output">Gathered output values</param>
            template <typename value_t>
            void allgather(const std::vector<value_t>& values, std::vector<value_t>& output) const;

            /// <summary>
            /// Allreduce a value per process
            /// </summary>
            /// <template name="value_t">Value type of the in- and output</template>
            /// <param name="value">Input value</param>
            /// <param name="output">Reduced output value</param>
            /// <param name="operation">Reduce operation</param>
            template <typename value_t>
            void allreduce(const value_t& value, value_t& output, MPI_Op operation) const;

            /// <summary>
            /// Allreduce in place a value per process
            /// </summary>
            /// <template name="value_t">Value type of the in- and output</template>
            /// <param name="value">In- and output value</param>
            /// <param name="operation">Reduce operation</param>
            template <typename value_t>
            void allreduce_inplace(value_t& value, MPI_Op operation) const;

            /// <summary>
            /// Barrier, waiting for all processes
            /// </summary>
            void barrier() const;

        private:
            /// <summary>
            /// Constructor
            /// </summary>
            mpi(MPI_Comm comm);

            /// Communicator
            MPI_Comm comm;

            /// Process information
            int rank;
            int num_processes;

            /// Was mpi already initialized before creating this object?
            bool already_initialized;
        };

        /// <summary>
        /// Singleton instance
        /// </summary>
        /// <return>mpi instance</return>
        static const mpi& get_instance(MPI_Comm comm = MPI_COMM_WORLD);
    }
}

#include "tpf_mpi.inl"