#include "tpf_loop.h"

#include <algorithm>
#include <iterator>
#include <type_traits>

namespace tpf
{
    namespace algorithm
    {
        namespace
        {
            template <bool check, typename function_t, typename arg_t>
            struct conditional_return
            {
                static typename std::result_of<function_t(arg_t)>::type call(function_t function, arg_t arg) { return function(arg); }
            };

            template <typename function_t, typename arg_t>
            struct conditional_return<false, function_t, arg_t>
            {
                static void* call(function_t function, arg_t arg) { function(arg); return nullptr; }
            };

            template <bool check, typename output_iterator_t>
            struct conditional_container
            {
                using type = typename output_iterator_t::container_type;
            };

            template <typename output_iterator_t>
            struct conditional_container<false, output_iterator_t>
            {
                using type = output_iterator_t;
            };

            template <bool check>
            struct conditional_inserter
            {
                template <typename container_t>
                static std::insert_iterator<container_t> make_inserter(container_t& container)
                {
                    return std::inserter(container, container.end());
                }
            };

            template <>
            struct conditional_inserter<false>
            {
                template <typename output_iterator_t>
                static output_iterator_t make_inserter(output_iterator_t& output_iterator)
                {
                    return output_iterator;
                }
            };

            template <bool check>
            struct conditional_copy
            {
                template <typename container_t, typename output_iterator_t>
                static void copy(const container_t& private_container, output_iterator_t out_begin)
                {
                    std::for_each(private_container.begin(), private_container.end(),
                        [&out_begin](const typename output_iterator_t::container_type::value_type& value) { out_begin++ = value; });
                }
            };

            template <>
            struct conditional_copy<false>
            {
                template <typename container_t, typename output_iterator_t>
                static void copy(const container_t&, output_iterator_t)
                {
                    return;
                }
            };
        }

        template <typename integral_t, int dimensions, typename function_t, typename output_iterator_t>
        inline void nested_loop(const Eigen::Matrix<integral_t, dimensions, 1>& start, const Eigen::Matrix<integral_t, dimensions, 1>& end,
            function_t function, output_iterator_t out_begin)
        {
            const auto max = (end - start).prod();

            for (integral_t index = 0; index < max; ++index)
            {
                // Calculate coordinates
                Eigen::Matrix<integral_t, dimensions, 1> coords;

                std::size_t index_rest = index;

                for (std::size_t dim = 0; dim < static_cast<std::size_t>(dimensions); ++dim)
                {
                    const std::size_t size = end[dim] - start[dim];

                    coords[dim] = start[dim] + index_rest % size;
                    index_rest /= size;
                }

                // Call function
                out_begin++ = conditional_return<!std::is_void<typename std::invoke_result<function_t, const Eigen::Matrix<integral_t, dimensions, 1>&>::type>::value,
                    function_t, const Eigen::Matrix<integral_t, dimensions, 1>&>::call(function, coords);
            }
        }

        template <typename integral_t, int dimensions, typename function_t, typename output_iterator_t>
        inline void parallel_nested_loop(const Eigen::Matrix<integral_t, dimensions, 1>& start, const Eigen::Matrix<integral_t, dimensions, 1>& end,
            function_t function, output_iterator_t out_begin)
        {
            using signed_integral_t = typename std::make_signed<integral_t>::type;

            const auto outer_start = static_cast<signed_integral_t>(start[dimensions - 1]);
            const auto outer_end = static_cast<signed_integral_t>(end[dimensions - 1]);

            #pragma omp parallel for schedule(dynamic)
            for (signed_integral_t outer_index = outer_start; outer_index < outer_end; ++outer_index)
            {
                // Create container given by the output iterator for temporarily storing private information
                constexpr bool dummy = std::is_same<output_iterator_t, dummy_output_iterator<typename std::invoke_result<function_t, Eigen::Matrix<integral_t, dimensions, 1>>::type>>::value;

                typename conditional_container<!dummy, output_iterator_t>::type private_container;

                // Create start and end for inner loops
                Eigen::Matrix<integral_t, dimensions, 1> inner_start, inner_end;

                inner_start << start.template head<dimensions - 1>(), outer_index;
                inner_end << end.template head<dimensions - 1>(), outer_index + 1;

                nested_loop(inner_start, inner_end, function, conditional_inserter<!dummy>::make_inserter(private_container));

                #pragma omp critical(write_particles)
                {
                    conditional_copy<!dummy>::copy(private_container, out_begin);
                }
            }
        }
    }
}
