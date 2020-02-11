#include "tpf_minmax.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

namespace tpf
{
    namespace utility
    {
        template <typename value_t>
        inline value_t min(value_t first, value_t second)
        {
            return std::min(first, second);
        }

        template <typename value_t, typename... other_t>
        inline value_t min(value_t first, value_t second, other_t... others)
        {
            return min(min(first, second), others...);
        }

        template <typename iterator_t, typename comparator>
        inline auto min_in_range(iterator_t first, iterator_t end, comparator comp) ->
            decltype(comp(std::declval<typename std::iterator_traits<iterator_t>::value_type>(), std::declval<typename std::iterator_traits<iterator_t>::value_type>()))
        {
            auto minimum = comp(*first, *first);

            for (auto current = ++first; current != end; ++current)
            {
                minimum = min(minimum, comp(*first, *current));
            }

            return minimum;
        }

        template <typename value_t>
        inline value_t max(value_t first, value_t second)
        {
            return std::max(first, second);
        }

        template <typename value_t, typename... other_t>
        inline value_t max(value_t first, value_t second, other_t... others)
        {
            return max(max(first, second), others...);
        }

        template <typename iterator_t, typename comparator>
        inline auto max_in_range(iterator_t first, iterator_t end, comparator comp) ->
            decltype(comp(std::declval<typename std::iterator_traits<iterator_t>::value_type>(), std::declval<typename std::iterator_traits<iterator_t>::value_type>()))
        {
            auto maximum = comp(*first, *first);

            for (auto current = ++first; current != end; ++current)
            {
                maximum = max(maximum, comp(*first, *current));
            }

            return maximum;
        }
    }
}
