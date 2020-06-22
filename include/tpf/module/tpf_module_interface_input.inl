#include "tpf_module_interface_input.h"

#include "tpf/log/tpf_log.h"

#include <stdexcept>

namespace tpf
{
    namespace modules
    {
        template <typename... Inputs>
        inline interface_input<Inputs...>::interface_input() : input_set(false) {}

        template <typename... Inputs>
        inline void interface_input<Inputs...>::set_input(Inputs... inputs)
        {
            set_algorithm_input(std::forward<Inputs>(inputs)...);

            this->input_set = true;
        }

        template <typename... Inputs>
        inline void interface_input<Inputs...>::sanity_check() const
        {
            if (!this->input_set)
            {
                throw std::runtime_error(__tpf_error_message("Input was not set properly."));
            }
        }
    }
}
