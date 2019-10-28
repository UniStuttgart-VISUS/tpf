#pragma once

#include <functional>
#include <type_traits>
#include <utility>

namespace tpf
{
    namespace utility
    {
        /// <summary>
        /// Minimum of two values
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <param name="first">First value</param>
        /// <param name="second">Second value</param>
        /// <returns>Minimum of both values</returns>
        template <typename value_t>
        value_t min(value_t first, value_t second);

        /// <summary>
        /// Minimum of multiple values
        /// </summary>
        /// <template name="value_t">Value type of first two values</template>
        /// <template name="other_t">Value type of subsequent values</template>
        /// <param name="first">First value</param>
        /// <param name="second">Second value</param>
        /// <param name="others">Other values</param>
        /// <returns>Minimum of values</returns>
        template <typename value_t, typename... other_t>
        value_t min(value_t first, value_t second, other_t... others);

        /// <summary>
        /// Minimum of multiple values
        /// </summary>
        /// <template name="iterator_t">Forward iterator type</template>
        /// <param name="begin">Iterator pointing at the first element in the range</param>
        /// <param name="end">Iterator pointing past the last element in the range</param>
        /// <returns>Minimum of values</returns>
        template <typename iterator_t, typename comparator = std::less<typename iterator_t::value_type>>
        auto min_in_range(iterator_t begin, iterator_t end, comparator comp = std::less<typename iterator_t::value_type>()) ->
            decltype(comp(std::declval<typename iterator_t::value_type>(), std::declval<typename iterator_t::value_type>()));

        /// <summary>
        /// Maximum of two values
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <param name="first">First value</param>
        /// <param name="second">Second value</param>
        /// <returns>Maximum of both values</returns>
        template <typename value_t>
        value_t max(value_t first, value_t second);

        /// <summary>
        /// Maximum of multiple values
        /// </summary>
        /// <template name="value_t">Value type of first two values</template>
        /// <template name="other_t">Value type of subsequent values</template>
        /// <param name="first">First value</param>
        /// <param name="second">Second value</param>
        /// <param name="others">Other values</param>
        /// <returns>Maximum of values</returns>
        template <typename value_t, typename... other_t>
        value_t max(value_t first, value_t second, other_t... others);

        /// <summary>
        /// Maximum of multiple values
        /// </summary>
        /// <template name="iterator_t">Forward iterator type</template>
        /// <param name="begin">Iterator pointing at the first element in the range</param>
        /// <param name="end">Iterator pointing past the last element in the range</param>
        /// <returns>Maximum of values</returns>
        template <typename iterator_t, typename comparator = std::less<typename iterator_t::value_type>>
        auto max_in_range(iterator_t begin, iterator_t end, comparator comp = std::greater<typename iterator_t::value_type>()) ->
            decltype(comp(std::declval<typename iterator_t::value_type>(), std::declval<typename iterator_t::value_type>()));
    }
}

#include "tpf_minmax.inl"
