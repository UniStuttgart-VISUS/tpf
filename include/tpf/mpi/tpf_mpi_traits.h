#pragma once

#ifdef __tpf_use_mpi
#include "mpi.h"
#endif

namespace tpf
{
    namespace mpi
    {
#ifdef __tpf_use_mpi
        /// Mapping of standard type to MPI type
        template <typename std_type> struct mpi_t { };
        template <> struct mpi_t<char> { static constexpr MPI_Datatype value = MPI_CHAR; typedef char type; };
        template <> struct mpi_t<short> { static constexpr MPI_Datatype value = MPI_SHORT; typedef short type; };
        template <> struct mpi_t<int> { static constexpr MPI_Datatype value = MPI_INT; typedef int type; };
        template <> struct mpi_t<long> { static constexpr MPI_Datatype value = MPI_LONG; typedef long type; };
        template <> struct mpi_t<float> { static constexpr MPI_Datatype value = MPI_FLOAT; typedef float type; };
        template <> struct mpi_t<double> { static constexpr MPI_Datatype value = MPI_DOUBLE; typedef double type; };
        template <> struct mpi_t<unsigned char> { static constexpr MPI_Datatype value = MPI_UNSIGNED_CHAR; typedef unsigned char type; };
        template <> struct mpi_t<unsigned short> { static constexpr MPI_Datatype value = MPI_UNSIGNED_SHORT; typedef unsigned short type; };
        template <> struct mpi_t<unsigned int> { static constexpr MPI_Datatype value = MPI_UNSIGNED; typedef unsigned int type; };
        template <> struct mpi_t<unsigned long> { static constexpr MPI_Datatype value = MPI_UNSIGNED_LONG; typedef unsigned long type; };
        template <> struct mpi_t<long double> { static constexpr MPI_Datatype value = MPI_LONG_DOUBLE; typedef long double type; };
        template <> struct mpi_t<long long> { static constexpr MPI_Datatype value = MPI_LONG_LONG_INT; typedef long long type; };

        /// Mapping of non-standard type to MPI type
        template <> struct mpi_t<unsigned long long> { static constexpr MPI_Datatype value = MPI_UNSIGNED_LONG; typedef unsigned long type; };
#endif

#ifndef __tpf_use_mpi
        /// Define communicator type and world communicator for non-MPI version
#define MPI_Comm void*
#define MPI_COMM_WORLD nullptr

        /// Define gather operation type and gather operations for non-MPI version
#define MPI_Op nullptr_t
#define MPI_OP_NULL nullptr
#define MPI_MAX     nullptr
#define MPI_MIN     nullptr
#define MPI_SUM     nullptr
#define MPI_PROD    nullptr
#define MPI_LAND    nullptr
#define MPI_BAND    nullptr
#define MPI_LOR     nullptr
#define MPI_BOR     nullptr
#define MPI_LXOR    nullptr
#define MPI_BXOR    nullptr
#define MPI_MINLOC  nullptr
#define MPI_MAXLOC  nullptr
#define MPI_REPLACE nullptr
#endif
    }
}