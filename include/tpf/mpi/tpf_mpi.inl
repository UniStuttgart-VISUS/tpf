#include "tpf_mpi.h"

#include "tpf_mpi_exceptions.h"
#include "tpf_mpi_traits.h"

#include "../log/tpf_log.h"

#include <algorithm>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        inline const mpi& mpi::get_instance(MPI_Comm comm)
        {
            struct make_shared_enabler : public mpi { public: make_shared_enabler(MPI_Comm comm) : mpi(comm) {} };

            static auto world_instance = std::static_pointer_cast<mpi>(std::make_shared<make_shared_enabler>(MPI_COMM_WORLD));
            static std::map<MPI_Comm, std::shared_ptr<mpi>> mpi_instances{ {world_instance->get_comm(), world_instance} };

            if (comm == MPI_COMM_WORLD)
            {
                return *world_instance;
            }
            else
            {
                if (mpi_instances.find(comm) == mpi_instances.end())
                {
                    mpi_instances[comm] = std::static_pointer_cast<mpi>(std::make_shared<make_shared_enabler>(comm));
                }

                return *mpi_instances.at(comm);
            }
        }

        inline mpi::mpi(MPI_Comm comm)
        {
            this->rank = 0;
            this->num_processes = 1;

#ifdef __tpf_use_mpi
            // Initialize mpi as needed
            int flag;
            MPI_Initialized(&flag);

            this->already_initialized = flag != 0;

            if (!this->already_initialized)
            {
                MPI_Init(nullptr, nullptr);
            }

            if (comm == MPI_COMM_WORLD)
            {
                MPI_Comm_dup(MPI_COMM_WORLD, &this->comm);
            }
            else
            {
                this->comm = comm;
            }

            // Get rank and number of processes
            MPI_Comm_rank(this->comm, &this->rank);
            MPI_Comm_size(this->comm, &this->num_processes);
#endif
        }

        inline mpi::~mpi()
        {
#ifdef __tpf_use_mpi
            if (!this->already_initialized)
            {
                MPI_Finalize();
            }
#endif
        }

        inline bool mpi::check_mpi_status() const
        {
            return get_num_processes() > 1;
        }

        inline int mpi::get_rank() const
        {
            return this->rank;
        }

        inline int mpi::get_num_processes() const
        {
            return this->num_processes;
        }

        inline MPI_Comm mpi::get_comm() const
        {
            return this->comm;
        }

        inline int mpi::communicate_error(const bool error) const
        {
            int errcode = error ? 0 : -1;

#ifdef __tpf_use_mpi
            if (check_mpi_status())
            {
                // Send error code
                const int send_errcode = error ? this->rank : -1;

                allreduce(send_errcode, errcode, MPI_MAX);
            }
#endif

            return errcode;
        }

        template <typename value_t>
        inline void mpi::broadcast(value_t& value, const int root) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Broadcast value
            auto buffer = static_cast<typename mpi_t<value_t>::type>(value);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            MPI_Bcast(&buffer, 1, mpi_t<value_t>::value, root, this->comm);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            value = static_cast<value_t>(buffer);
#endif
        }

        template <typename value_t>
        inline void mpi::broadcast(const value_t& value, const int root) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
#ifdef __tpf_sanity_checks
            // Check that it is sending and not receiving
            if (root != this->rank)
            {
                throw mpi_exception(__tpf_error_message("Broadcasting a constant object only allowed for sender"));
            }
#endif

            broadcast(const_cast<value_t&>(value), root);
#endif
        }

        template <typename value_t>
        inline void mpi::broadcast(std::vector<value_t>& values, const int root) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Broadcast array size
            std::size_t buffer_size = values.size();

            broadcast(buffer_size, root);

            // Broadcast data
            if (std::is_same<typename mpi_t<value_t>::type, value_t>::value)
            {
                values.resize(buffer_size);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Bcast(values.data(), static_cast<int>(buffer_size), mpi_t<value_t>::value, root, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
            else
            {
                std::vector<typename mpi_t<value_t>::type> buffer(buffer_size);

                if (root == this->rank)
                {
                    std::transform(values.begin(), values.end(), buffer.begin(),
                        [](const value_t& value) { return static_cast<typename mpi_t<value_t>::type>(value); });
                }

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Bcast(buffer.data(), static_cast<int>(buffer_size), mpi_t<value_t>::value, root, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                if (root != this->rank)
                {
                    values.resize(buffer_size);

                    std::transform(buffer.begin(), buffer.end(), values.begin(),
                        [](const typename mpi_t<value_t>::type& value) { return static_cast<value_t>(value); });
                }
            }
#endif
        }

        template <typename value_t>
        inline void mpi::broadcast(const std::vector<value_t>& values, const int root) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
#ifdef __tpf_sanity_checks
            // Check that it is sending and not receiving
            if (root != this->rank)
            {
                throw mpi_exception(__tpf_error_message("Broadcasting a constant object only allowed for sender"));
            }
#endif

            broadcast(const_cast<std::vector<value_t>&>(values), root);
#endif
        }

        template <typename value_t>
        inline void mpi::send(const value_t& value, const int destination) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Send value
            auto buffer = static_cast<typename mpi_t<value_t>::type>(value);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            MPI_Send(&buffer, 1, mpi_t<value_t>::value, destination, 0, this->comm);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif
#endif
        }

        template <typename value_t>
        inline void mpi::send(const std::vector<value_t>& values, const int destination) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Send array size
            send(values.size(), destination);

            // Send data
            if (std::is_same<typename mpi_t<value_t>::type, value_t>::value)
            {
#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Send(const_cast<value_t*>(values.data()), static_cast<int>(values.size()), mpi_t<value_t>::value, destination, 0, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
            else
            {
                std::vector<typename mpi_t<value_t>::type> buffer(values.size());

                std::transform(values.begin(), values.end(), buffer.begin(),
                    [](const value_t& value) { return static_cast<typename mpi_t<value_t>::type>(value); });

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Send(buffer.data(), static_cast<int>(values.size()), mpi_t<value_t>::value, destination, 0, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
#endif
        }

        template <typename value_t>
        inline void mpi::receive(value_t& value, const int source) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Receive value
            typename mpi_t<value_t>::type buffer;

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            MPI_Recv(&buffer, 1, mpi_t<value_t>::value, source, 0, this->comm, MPI_STATUS_IGNORE);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            value = static_cast<value_t>(buffer);
#endif
        }

        template <typename value_t>
        inline void mpi::receive(std::vector<value_t>& values, const int source) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // Receive array size
            std::size_t buffer_size;

            receive(buffer_size, source);

            // Receive data
            values.resize(buffer_size);

            if (std::is_same<typename mpi_t<value_t>::type, value_t>::value)
            {
#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Recv(values.data(), static_cast<int>(values.size()), mpi_t<value_t>::value, source, 0, this->comm, MPI_STATUS_IGNORE);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
            else
            {
                std::vector<typename mpi_t<value_t>::type> buffer(values.size());

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Recv(buffer.data(), static_cast<int>(values.size()), mpi_t<value_t>::value, source, 0, this->comm, MPI_STATUS_IGNORE);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                std::transform(buffer.begin(), buffer.end(), values.begin(),
                    [](const typename mpi_t<value_t>::type& value) { return static_cast<value_t>(value); });
            }
#endif
        }

        template <typename value_t>
        inline void mpi::allgather(const value_t& value, std::vector<value_t>& output) const
        {
            // Return input if MPI is not used
            if (!check_mpi_status())
            {
                output.resize(1);
                output[0] = value;
                return;
            }

#ifdef __tpf_use_mpi
            // Allgather value
            if (std::is_same<typename mpi_t<value_t>::type, value_t>::value)
            {
                output.resize(this->num_processes);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Allgather(const_cast<value_t*>(&value), 1, mpi_t<value_t>::value, output.data(), 1, mpi_t<value_t>::value, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
            else
            {
                auto send_buffer = static_cast<typename mpi_t<value_t>::type>(value);
                std::vector<typename mpi_t<value_t>::type> receive_buffer(this->num_processes);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Allgather(&send_buffer, 1, mpi_t<value_t>::value, receive_buffer.data(), 1, mpi_t<value_t>::value, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                std::transform(receive_buffer.begin(), receive_buffer.end(), output.begin(),
                    [](const typename mpi_t<value_t>::type& value) { return static_cast<value_t>(value); });
            }
#endif
        }

        template <typename value_t>
        inline void mpi::allgather(const std::vector<value_t>& values, std::vector<value_t>& output) const
        {
            // Return input if MPI is not used
            if (!check_mpi_status())
            {
                output = values;
                return;
            }

#ifdef __tpf_use_mpi
            // Allgather value
            if (std::is_same<typename mpi_t<value_t>::type, value_t>::value)
            {
                output.resize(this->num_processes * values.size());

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Allgather(const_cast<value_t*>(values.data()), static_cast<int>(values.size()), mpi_t<value_t>::value,
                    output.data(), static_cast<int>(values.size()), mpi_t<value_t>::value, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif
            }
            else
            {
                std::vector<typename mpi_t<value_t>::type> send_buffer(values.size());
                std::vector<typename mpi_t<value_t>::type> receive_buffer(this->num_processes * values.size());

                std::transform(values.begin(), values.end(), send_buffer.begin(),
                    [](const value_t& value) { return static_cast<typename mpi_t<value_t>::type>(value); });

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                MPI_Allgather(&send_buffer, static_cast<int>(send_buffer.size()), mpi_t<value_t>::value, receive_buffer.data(),
                    static_cast<int>(send_buffer.size()), mpi_t<value_t>::value, this->comm);

#ifdef __tpf_debug
                MPI_Barrier(this->comm);
#endif

                std::transform(receive_buffer.begin(), receive_buffer.end(), output.begin(),
                    [](const typename mpi_t<value_t>::type& value) { return static_cast<value_t>(value); });
            }
#endif
        }

        template <typename value_t>
        inline void mpi::allreduce(const value_t& value, value_t& output, MPI_Op operation) const
        {
            // Return input if MPI is not used
            if (!check_mpi_status())
            {
                output = value;
                return;
            }

#ifdef __tpf_use_mpi
            // All reduce
            auto send_buffer = static_cast<typename mpi_t<value_t>::type>(value);
            typename mpi_t<value_t>::type receive_buffer;

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            MPI_Allreduce(&send_buffer, &receive_buffer, 1, mpi_t<value_t>::value, operation, this->comm);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            output = static_cast<value_t>(receive_buffer);
#endif
        }

        template <typename value_t>
        inline void mpi::allreduce_inplace(value_t& value, MPI_Op operation) const
        {
            // Return if MPI is not used
            if (!check_mpi_status())
            {
                return;
            }

#ifdef __tpf_use_mpi
            // All reduce
            auto buffer = static_cast<typename mpi_t<value_t>::type>(value);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            MPI_Allreduce(MPI_IN_PLACE, &buffer, 1, mpi_t<value_t>::value, operation, this->comm);

#ifdef __tpf_debug
            MPI_Barrier(this->comm);
#endif

            value = static_cast<value_t>(buffer);
#endif
        }

        inline void mpi::barrier() const
        {
#ifdef __tpf_use_mpi
            MPI_Barrier(this->comm);
#endif
        }

        inline const mpi& get_instance(MPI_Comm comm)
        {
            return mpi::get_instance(comm);
        }
    }
}
