#include "tpf_math_defines.h"

#include <cmath>

namespace tpf
{
    namespace math
    {
        template <typename floatp_t>
        const floatp_t default_epsilon<floatp_t>::value = static_cast<floatp_t>(0.0L);

        template <typename floatp_t>
        inline constexpr typename std::enable_if<std::is_floating_point<floatp_t>::value, bool>::type equals(const floatp_t first, const floatp_t second, const floatp_t epsilon)
        {
            return std::abs(first - second) <= epsilon;
        }
    }
}