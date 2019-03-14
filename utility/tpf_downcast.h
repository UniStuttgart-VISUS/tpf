#pragma once

#include <limits>
#include <type_traits>

namespace tpf
{
    namespace utility
    {
        namespace
        {
            template <typename new_t, typename original_t>
            struct narrowing
            {
                static constexpr bool value = 
                    static_cast<original_t>(std::numeric_limits<new_t>::max()) < std::numeric_limits<original_t>::max() ||
                    static_cast<original_t>(std::numeric_limits<new_t>::lowest()) > std::numeric_limits<original_t>::lowest();
            };

            template <typename new_t, typename original_t>
            struct widening
            {
                static constexpr bool value = !(narrowing<new_t, original_t>::value || std::is_same<new_t, original_t>::value);
            };
        }

        /// <summary>
        /// Performs a cast with sanity checks for safe downcasting into (possibly) narrowing types
        /// </summary>
        /// <template name="new_t">New type</template>
        /// <template name="original_t">Original type</template>
        /// <param name="original">Original value</param>
        /// <returns>Converted value</returns>
        template <typename new_t, typename original_t>
        new_t safe_downcast(original_t original, typename std::enable_if<narrowing<new_t, original_t>::value>::type* = nullptr);

        /// <summary>
        /// Returns the original value, as it is of the same type
        /// </summary>
        /// <template name="new_t">New type</template>
        /// <template name="original_t">Original type</template>
        /// <param name="original">Original value</param>
        /// <returns>Converted value</returns>
        template <typename new_t, typename original_t>
        new_t safe_downcast(const original_t& original, typename std::enable_if<std::is_same<new_t, original_t>::value>::type* = nullptr);

        /// <summary>
        /// Casts and returns the original value, as it is not a narrowing cast
        /// </summary>
        /// <template name="new_t">New type</template>
        /// <template name="original_t">Original type</template>
        /// <param name="original">Original value</param>
        /// <returns>Converted value</returns>
        template <typename new_t, typename original_t>
        new_t safe_downcast(const original_t& original, typename std::enable_if<widening<new_t, original_t>::value>::type* = nullptr);
    }
}

#include "tpf_downcast.inl"