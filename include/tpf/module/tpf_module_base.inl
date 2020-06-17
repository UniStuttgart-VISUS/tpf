#include "tpf_module_base.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"
#include "tpf/mpi/tpf_mpi_exceptions.h"

#include <exception>
#include <stdexcept>

namespace tpf
{
    namespace modules
    {
        template <typename Parameters, typename Callbacks>
        inline module_base<Parameters, Callbacks>::module_base() : parameters_set(false), callbacks_set(false)
        {
        }

        template <typename Parameters, typename Callbacks>
        inline void module_base<Parameters, Callbacks>::run()
        {
            // Perform sanity checks on parameters and callbacks
            if (!sanity_check())
            {
                tpf::mpi::get_instance().communicate_error(true);

                if (!this->parameters_set && !this->callbacks_set)
                {
                    throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed. ",
                        "Parameters and callbacks were not set properly."));
                }
                else if (!this->parameters_set)
                {
                    throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed. ",
                        "Parameters were not set properly."));
                }
                else if (!this->callbacks_set)
                {
                    throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed. ",
                        "Callbacks were not set properly."));
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed. ",
                        "Sanity checks on parameters and callbacks revealed an error."));
                }
            }

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
            catch (const std::runtime_error& ex)
            {
                tpf::mpi::get_instance().communicate_error(true);

                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Module run of '", get_name(), "' failed."));
            }
            catch (const std::exception&)
            {
                tpf::mpi::get_instance().communicate_error(true);

                throw std::runtime_error(__tpf_error_message("Module run of '", get_name(), "' failed due to an unexpected error. ",
                    "This was probably caused by a programming mistake or unsufficient system resources. ",
                    "Try enabling TPF_OPTION_SANITY_CHECKS and TPF_OPTION_RANGE_CHECKS for higher security and better error handling."));
            }
        }

        template <typename Parameters, typename Callbacks>
        template <typename param_t>
        inline void module_base<Parameters, Callbacks>::set_parameters(typename std::enable_if<!std::is_same<param_t, void>::value, parameters_t>::type parameters)
        {
            set_algorithm_parameters(std::forward<parameters_t>(parameters));

            this->parameters_set = true;
        }

        template <typename Parameters, typename Callbacks>
        template <typename param_t>
        inline void module_base<Parameters, Callbacks>::set_parameters(typename std::enable_if<std::is_same<param_t, void>::value>::type)
        {
            this->parameters_set = true;
        }

        template <typename Parameters, typename Callbacks>
        template <typename cb_t>
        inline void module_base<Parameters, Callbacks>::set_callbacks(typename std::enable_if<!std::is_same<cb_t, void>::value, callbacks_t>::type callbacks)
        {
            set_algorithm_callbacks(std::forward<callbacks_t>(callbacks));

            this->callbacks_set = true;
        }

        template <typename Parameters, typename Callbacks>
        template <typename cb_t>
        inline void module_base<Parameters, Callbacks>::set_callbacks(typename std::enable_if<std::is_same<cb_t, void>::value>::type)
        {
            this->callbacks_set = true;
        }

        template <typename Parameters, typename Callbacks>
        inline void module_base<Parameters, Callbacks>::set_algorithm_parameters(parameters_t)
        {
        }

        template <typename Parameters, typename Callbacks>
        inline void module_base<Parameters, Callbacks>::set_algorithm_callbacks(callbacks_t)
        {
        }

        template <typename Parameters, typename Callbacks>
        inline bool module_base<Parameters, Callbacks>::sanity_check() const
        {
            return this->parameters_set && this->callbacks_set;
        }
    }
}