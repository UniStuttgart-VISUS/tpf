#pragma once

#include <exception>
#include <stdexcept>

namespace tpf
{
    namespace mpi
    {
        /// <summary>
        /// Exception for aborting all other mpi tasks cleanly
        /// </summary>
        class mpi_abort : public std::exception
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="rank">Rank of the failed process</param>
            mpi_abort(int rank);

            /// <summary>
            /// Return the rank of the failed process
            /// </summary>
            /// <returns>Rank of the failed process</returns>
            int get_rank() const;

        private:
            /// Rank of the failed process
            int rank;
        };

        /// <summary>
        /// Exception for handling mpi errors
        /// </summary>
        class mpi_exception : public std::runtime_error
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="message">Error message</param>
            mpi_exception(const std::string& message);
        };
    }
}

#include "tpf_mpi_exceptions.inl"