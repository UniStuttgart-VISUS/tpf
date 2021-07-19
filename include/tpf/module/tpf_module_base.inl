#include "tpf_module_base.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"
#include "tpf/mpi/tpf_mpi_exceptions.h"

#include <exception>
#include <functional>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace tpf
{
    namespace modules
    {
        template <typename... Interfaces>
        inline module_base<Interfaces...>::module_base()
        {
        }

        template <typename... Interfaces>
        inline void module_base<Interfaces...>::run()
        {
            // Perform sanity checks on input, output, parameters and callbacks and throw
            // exception if failed
            sanity_checks();

            // Run the algorithm
            try
            {
                run_algorithm();

                if (int errcode = mpi::get_instance().communicate_error() != -1)
                {
                    throw mpi::mpi_abort(errcode);
                }
            }
            catch (const mpi::mpi_abort& ex)
            {
                throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' aborted. Exception thrown on rank ", ex.get_rank()));
            }
            catch (const mpi::mpi_exception& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Module run of '", get_name(), "' failed due to MPI error."));
            }
            catch (const std::exception& ex)
            {
                tpf::mpi::get_instance().communicate_error(true);

                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Module run of '", get_name(), "' failed: ", ex.what()));
            }
            catch (...)
            {
                tpf::mpi::get_instance().communicate_error(true);

                throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed due to an unexpected error. ",
                    "This was probably caused by a programming mistake or unsufficient system resources. ",
                    "Try enabling TPF_OPTION_SANITY_CHECKS and TPF_OPTION_RANGE_CHECKS for higher security and better error handling."));
            }
        }

        template<typename... Interfaces>
        inline void module_base<Interfaces...>::sanity_checks() const
        {
            int dummy[] = { (Interfaces::sanity_check(), 0)... };
        }

        template <typename T>
        inline T get_or_default(std::optional<T> parameter, T default_value)
        {
            if (parameter)
            {
                return *parameter;
            }

            return default_value;
        }

        template <typename T>
        inline T* get_or_default(std::optional<std::reference_wrapper<T>> parameter, T* default_value)
        {
            if (parameter)
            {
                return &parameter->get();
            }

            return default_value;
        }

        template <typename T>
        inline std::optional<std::reference_wrapper<T>> set_or_default(T& parameter, bool set)
        {
            if (set)
            {
                return parameter;
            }

            return std::nullopt;
        }

        template <typename T>
        inline std::optional<std::reference_wrapper<T>> set_or_default(T* parameter)
        {
            if (parameter != nullptr)
            {
                return *parameter;
            }

            return std::nullopt;
        }
    }
}
