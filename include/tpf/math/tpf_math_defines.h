#pragma once

#include <type_traits>

namespace tpf
{
    namespace math
    {
        /// pi
        template <typename floatp_t> static const floatp_t pi = static_cast<floatp_t>(std::acos(-1.0L));

        /// default epsilon for errors in functions
        template <typename floatp_t>
        struct default_epsilon
        {
            static const floatp_t value;
        };

        /// <summary>
        /// Comparison function for floating point types with error margin
        /// </summary>
        /// <template name="floatp_t">Floating point type for values</template>
        /// <param name="first">First floating point value</param>
        /// <param name="second">Second floating point value</param>
        /// <param name="epsilon">Error margin</param>
        /// <returns>True, if values are within a certain margin of each other; false otherwise</returns>
        template <typename floatp_t>
        constexpr typename std::enable_if<std::is_floating_point<floatp_t>::value, bool>::type equals(floatp_t first, floatp_t second, floatp_t epsilon = default_epsilon<floatp_t>::value);
    }
}

#include "tpf_math_defines.inl"