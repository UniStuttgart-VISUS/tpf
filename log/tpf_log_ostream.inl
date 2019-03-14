#include "tpf_log_ostream.h"

#ifdef __tpf_use_mpi
#include "mpi.h"
#endif

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace tpf
{
    namespace log
    {
        template <typename error_stream_t, typename warning_stream_t, typename info_stream_t>
        inline log_ostream<error_stream_t, warning_stream_t, info_stream_t>::log_ostream(error_stream_t& err_stream, warning_stream_t& warn_stream, info_stream_t& info_stream)
            : log_base(), err_stream(err_stream), warn_stream(warn_stream), info_stream(info_stream), rank(0), num_processes(1)
#ifdef __tpf_use_mpi
            , world_comm(MPI_COMM_WORLD)
#endif
        { }

        template <typename error_stream_t, typename warning_stream_t, typename info_stream_t>
        inline void log_ostream<error_stream_t, warning_stream_t, info_stream_t>::write_error_message(const std::string& message, const bool root_only) const
        {
            write_message(this->err_stream, message, root_only);
        }

        template <typename error_stream_t, typename warning_stream_t, typename info_stream_t>
        inline void log_ostream<error_stream_t, warning_stream_t, info_stream_t>::write_warning_message(const std::string& message, const bool root_only) const
        {
            write_message(this->warn_stream, message, root_only);
        }

        template <typename error_stream_t, typename warning_stream_t, typename info_stream_t>
        inline void log_ostream<error_stream_t, warning_stream_t, info_stream_t>::write_info_message(const std::string& message, const bool root_only) const
        {
            write_message(this->info_stream, message, root_only);
        }

        template <typename error_stream_t, typename warning_stream_t, typename info_stream_t>
        template <typename stream_t>
        inline void log_ostream<error_stream_t, warning_stream_t, info_stream_t>::write_message(stream_t& stream, const std::string& message, const bool root_only) const
        {
#ifdef __tpf_use_mpi
            int flag;
            MPI_Initialized(&flag);

            if (flag == 0)
            {
                MPI_Init(nullptr, nullptr);
            }

            if (this->world_comm == MPI_COMM_WORLD)
            {
                MPI_Comm_dup(MPI_COMM_WORLD, &this->world_comm);

                MPI_Comm_rank(this->world_comm, &this->rank);
                MPI_Comm_size(this->world_comm, &this->num_processes);
            }

            const auto number_width = static_cast<int>(std::ceil(std::log10(this->num_processes)));
#endif

            if (root_only)
            {
                if (this->rank == 0)
                {
                    std::stringstream ss;
#ifdef __tpf_use_mpi
                    if (this->num_processes > 1)
                    {
                        ss << std::setfill('0') << std::setw(number_width) << this->rank;
                        ss << ": ";
                    }
#endif
                    ss << message << std::endl;

                    stream << ss.str();
                    stream.flush();
                }
            }
            else
            {
                for (int i = 0; i < this->num_processes; ++i)
                {
#ifdef __tpf_use_mpi
                    MPI_Barrier(this->world_comm);
#endif

                    if (i == this->rank)
                    {
                        std::stringstream ss;
#ifdef __tpf_use_mpi
                        if (this->num_processes > 1)
                        {
                            ss << std::setfill('0') << std::setw(number_width) << this->rank;
                            ss << ": ";
                        }
#endif
                        ss << message << std::endl;

                        stream << ss.str();
                        stream.flush();
                    }
                }

#ifdef __tpf_use_mpi
                MPI_Barrier(this->world_comm);
#endif
            }
        }
    }
}