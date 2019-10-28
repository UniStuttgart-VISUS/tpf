#include "tpf_not_implemented_exception.h"

#include <stdexcept>

namespace tpf
{
    namespace exception
    {
        inline not_implemented_exception::not_implemented_exception(std::string message)
            : std::logic_error(message.empty() ? "Function not yet implemented" : std::move(message)) { };
    }
}