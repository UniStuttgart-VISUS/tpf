#include "tpf_downcast.h"

#include "../log/tpf_log.h"

#include <limits>
#include <stdexcept>
#include <type_traits>

namespace tpf
{
    namespace utility
    {
        template <typename new_t, typename original_t>
        inline new_t safe_downcast(original_t original, typename std::enable_if<narrowing<new_t, original_t>::value>::type*)
        {
            if (!std::is_signed<new_t>::value && std::is_signed<original_t>::value && original < 0)
            {
                throw std::runtime_error(__tpf_error_message("Illegal conversion from negative value to unsigned type."));
            }
            else if (original > static_cast<original_t>(std::numeric_limits<new_t>::max()))
            {
                throw std::runtime_error(__tpf_error_message("Illegal narrowing conversion from large positive value."));
            }
            else if (std::is_signed<new_t>::value && std::is_signed<original_t>::value && original < static_cast<original_t>(std::numeric_limits<new_t>::lowest()))
            {
                throw std::runtime_error(__tpf_error_message("Illegal narrowing conversion from large negative value."));
            }

            return static_cast<new_t>(original);
        }

        template <typename new_t, typename original_t>
        inline new_t safe_downcast(const original_t& original, typename std::enable_if<std::is_same<new_t, original_t>::value>::type*)
        {
            return original;
        }

        template <typename new_t, typename original_t>
        inline new_t safe_downcast(const original_t& original, typename std::enable_if<widening<new_t, original_t>::value>::type*)
        {
            return static_cast<new_t>(original);
        }
    }
}