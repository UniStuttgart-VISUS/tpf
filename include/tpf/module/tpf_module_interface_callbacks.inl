#include "tpf_module_interface_callbacks.h"

#include "tpf/log/tpf_log.h"

#include <stdexcept>

namespace tpf
{
    namespace modules
    {
        template <typename... Callbacks>
        inline interface_callbacks<Callbacks...>::interface_callbacks() : callbacks_set(false) {}

        template <typename... Callbacks>
        inline void interface_callbacks<Callbacks...>::set_callbacks(Callbacks... callbacks)
        {
            set_algorithm_callbacks(std::forward<Callbacks>(callbacks)...);

            this->callbacks_set = true;
        }

        template <typename... Callbacks>
        inline void interface_callbacks<Callbacks...>::sanity_check() const
        {
            if (!this->callbacks_set)
            {
                throw std::runtime_error(__tpf_error_message("Callbacks were not set properly."));
            }
        }
    }
}
