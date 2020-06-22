#include "tpf_module_interface_output.h"

#include "tpf/log/tpf_log.h"

#include <stdexcept>

namespace tpf
{
    namespace modules
    {
        template <typename... Outputs>
        inline interface_output<Outputs...>::interface_output() : output_set(false) {}

        template <typename... Outputs>
        inline void interface_output<Outputs...>::set_output(Outputs... outputs)
        {
            set_algorithm_output(std::forward<Outputs>(outputs)...);

            this->output_set = true;
        }

        template <typename... Outputs>
        inline void interface_output<Outputs...>::sanity_check() const
        {
            if (!this->output_set)
            {
                throw std::runtime_error(__tpf_error_message("Output was not set properly."));
            }
        }
    }
}
