#include "tpf_mpi_exceptions.h"

#include <exception>
#include <stdexcept>

inline tpf::mpi::mpi_abort::mpi_abort(const int rank) : std::exception(), rank(rank) {};

inline int tpf::mpi::mpi_abort::get_rank() const
{
    return this->rank;
}

inline tpf::mpi::mpi_exception::mpi_exception(const std::string& message) : std::runtime_error(message) {};