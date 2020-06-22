#include "tpf_module_interface_parameters.h"

#include "tpf/log/tpf_log.h"

#include <stdexcept>

namespace tpf
{
    namespace modules
    {
        template <typename... Parameters>
        inline interface_parameters<Parameters...>::interface_parameters() : parameters_set(false) {}

        template <typename... Parameters>
        inline void interface_parameters<Parameters...>::set_parameters(Parameters... parameters)
        {
            set_algorithm_parameters(std::forward<Parameters>(parameters)...);

            this->parameters_set = true;
        }

        template <typename... Parameters>
        inline void interface_parameters<Parameters...>::sanity_check() const
        {
            if (!this->parameters_set)
            {
                throw std::runtime_error(__tpf_error_message("Parameters were not set properly."));
            }
        }
    }
}
