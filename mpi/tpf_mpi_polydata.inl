#include "tpf_mpi_polydata.h"

#include "tpf_mpi.h"
#include "tpf_mpi_exceptions.h"
#include "tpf_mpi_traits.h"

#include "../data/tpf_array.h"
#include "../data/tpf_array_base.h"
#include "../data/tpf_polydata.h"
#include "../data/tpf_data_information.h"

#include "../geometry/tpf_factory.h"
#include "../geometry/tpf_geometric_object.h"

#include "../log/tpf_log.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace mpi
    {
        template <typename point_t, typename kernel_t, typename... data_information_t>
        inline data::polydata<point_t> mpi_polydata::all_gather(const data::polydata<point_t>& source, const MPI_Comm comm, const data_information_t&... data_infos)
        {
#ifdef __tpf_use_mpi
            const auto& mpi = get_instance(comm);

            // Sanity checks
            if (!mpi.check_mpi_status())
            {
                return source;
            }

            // Test for exceptions before continuing
            if (int errcode = mpi.communicate_error() != -1)
            {
                throw mpi_abort(errcode);
            }

            // In turn send and receive data
            data::polydata<point_t> target;

            for (int proc = 0; proc < mpi.get_num_processes(); ++proc)
            {
                if (proc == mpi.get_rank())
                {
                    // Broadcast
                    send_polydata(source, proc, true, MPI_COMM_WORLD, data_infos...);
                    target.merge(source);
                }
                else
                {
                    // Receive broadcast
                    auto received = receive_polydata<point_t, kernel_t>(proc, true, MPI_COMM_WORLD, data_infos...);
                    target.merge(received);
                }
            }

            return target;
#else
            return source;
#endif
        }

#ifdef __tpf_use_mpi
        template <typename point_t>
        inline void mpi_polydata::send_polydata(const data::polydata<point_t>& source, const int receiver, const bool broadcast, const MPI_Comm comm)
        {
            // Send geometry only
            send_geometry(source->get_geometry(), receiver, broadcast, comm);
        }

        template <typename point_t, typename... data_information_t>
        inline void mpi_polydata::send_polydata(const data::polydata<point_t>& source, const int receiver, const bool broadcast, const MPI_Comm comm, const data_information_t&... data_infos)
        {
            // First send geometry
            send_geometry(source.get_geometry(), receiver, broadcast, comm);

            // Then send attached arrays
            send_polydata_arrays(source, receiver, broadcast, comm, data_infos...);
        }

        template <typename point_t>
        inline void mpi_polydata::send_geometry(const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& source, const int receiver, const bool broadcast, MPI_Comm comm)
        {
            static_assert(std::is_same<point_t, typename mpi_t<point_t>::type>::value, "Value type must be mpi_polydata compatible");

            const auto& mpi = get_instance(comm);

            // Send number of objects
            if (broadcast)
            {
                mpi.broadcast(source.size(), mpi.get_rank());
            }
            else
            {
                mpi.broadcast(source.size(), receiver);
            }

            // For each object
            for (const auto object : source)
            {
                // Serialize object
                auto serialized = object->serialize();

                // Send object
                if (broadcast)
                {
                    mpi.broadcast(serialized, mpi.get_rank());
                }
                else
                {
                    mpi.send(serialized, receiver);
                }
            }
        }

        template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_information_t>
        inline void mpi_polydata::send_polydata_arrays(const data::polydata<point_t>& source, const int receiver, const bool broadcast, const MPI_Comm comm,
            const data::data_information<value_t, rows, columns>& first_data_info, const data_information_t&... more_data_infos)
        {
            send_polydata_arrays(source, receiver, broadcast, comm, first_data_info);
            send_polydata_arrays(source, receiver, broadcast, comm, more_data_infos...);
        }

        template <typename point_t, typename value_t, std::size_t rows, std::size_t columns>
        inline void mpi_polydata::send_polydata_arrays(const data::polydata<point_t>& source, const int receiver, const bool broadcast, const MPI_Comm comm,
            const data::data_information<value_t, rows, columns>& data_info)
        {
            const auto& mpi = get_instance(comm);

            // Send array name
            const std::vector<std::string::value_type> name(data_info.name.begin(), data_info.name.end());

            if (broadcast)
            {
                mpi.broadcast(name, mpi.get_rank());
            }
            else
            {
                mpi.send(name, receiver);
            }

            // Send array
            if (data_info.topology == data::topology_t::POINT_DATA)
            {
                mpi_base::send_array(source.template get_point_data_as<value_t, rows, columns>(data_info.name)->get_data(), receiver, broadcast, comm);
            }
            else if (data_info.topology == data::topology_t::CELL_DATA)
            {
                mpi_base::send_array(source.template get_cell_data_as<value_t, rows, columns>(data_info.name)->get_data(), receiver, broadcast, comm);
            }
            else if (data_info.topology == data::topology_t::OBJECT_DATA)
            {
                mpi_base::send_array(source.template get_object_data_as<value_t, rows, columns>(data_info.name)->get_data(), receiver, broadcast, comm);
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }
        }

        template <typename point_t, typename kernel_t>
        inline data::polydata<point_t> mpi_polydata::receive_polydata(const int sender, const bool broadcast, const MPI_Comm comm)
        {
            // Receive only geometry
            data::polydata<point_t> received(receive_geometry<point_t, kernel_t>(sender, broadcast, comm));

            return received;
        }

        template <typename point_t, typename kernel_t, typename... data_information_t>
        inline data::polydata<point_t> mpi_polydata::receive_polydata(const int sender, const bool broadcast, const MPI_Comm comm, const data_information_t&... data_infos)
        {
            // First receive geometry
            data::polydata<point_t> received(receive_geometry<point_t, kernel_t>(sender, broadcast, comm));

            // Then receive attached arrays
            receive_polydata_arrays(received, sender, broadcast, comm, data_infos...);

            return received;
        }

        template <typename point_t, typename kernel_t>
        inline std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> mpi_polydata::receive_geometry(const int sender, const bool broadcast, MPI_Comm comm)
        {
            static_assert(std::is_same<point_t, typename mpi_t<point_t>::type>::value, "Value type must be mpi_polydata compatible");

            const auto& mpi = get_instance(comm);

            // Receive number of objects
            typename mpi_t<std::size_t>::type num_objects;

            if (broadcast)
            {
                mpi.broadcast(num_objects, sender);
            }
            else
            {
                mpi.receive(num_objects, sender);
            }

            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> objects;
            objects.reserve(num_objects);

            // For each object
            for (typename mpi_t<std::size_t>::type i = 0; i < num_objects; ++i)
            {
                // Receive object
                std::vector<char> data;

                if (broadcast)
                {
                    mpi.broadcast(data, sender);
                }
                else
                {
                    mpi.receive(data, sender);
                }

                // Deserialize object
                objects.push_back(geometry::template deserialize<point_t, kernel_t>(data));
            }

            return objects;
        }

        template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_information_t>
        inline void mpi_polydata::receive_polydata_arrays(data::polydata<point_t>& target, const int sender, const bool broadcast, const MPI_Comm comm,
            const data::data_information<value_t, rows, columns>& first_data_info, const data_information_t&... more_data_infos)
        {
            receive_polydata_arrays(target, sender, broadcast, comm, first_data_info);
            receive_polydata_arrays(target, sender, broadcast, comm, more_data_infos...);
        }

        template <typename point_t, typename value_t, std::size_t rows, std::size_t columns>
        inline void mpi_polydata::receive_polydata_arrays(data::polydata<point_t>& target, const int sender, const bool broadcast, const MPI_Comm comm,
            const data::data_information<value_t, rows, columns>& data_info)
        {
            const auto& mpi = get_instance(comm);

            // Receive array name
            std::vector<typename std::string::value_type> buffer;

            if (broadcast)
            {
                mpi.broadcast(buffer, sender);
            }
            else
            {
                mpi.receive(buffer, sender);
            }

            std::string name(buffer.begin(), buffer.end());

            // Receive array
            auto data = mpi_base::receive_array<float>(sender, broadcast, comm);

            target.add(std::dynamic_pointer_cast<data::array_base>(std::make_shared<data::array<value_t, rows, columns>>(data_info.name, std::move(data))), data_info.topology);
        }
#endif
    }
}
