#pragma once

#include "Eigen/Dense"

#include <type_traits>

namespace tpf
{
    namespace algorithm
    {
        /// <summary>
        /// Dummy output iterator used as default with no effect
        /// </summary>
        template <typename value_t>
        struct dummy_output_iterator
        {
        public:
            dummy_output_iterator() = default;
            dummy_output_iterator(const dummy_output_iterator&) = default;

            dummy_output_iterator& operator=(const dummy_output_iterator&) = default;

            template <typename local_value_t = value_t>
            dummy_output_iterator& operator=(void*) { return *this; }

            dummy_output_iterator& operator++() { return *this; }
            dummy_output_iterator operator++(int) { return *this; }
        };

        /// <summary>
        /// Dynamic nested loop for all permutations in [{start_x, start_y, ...}, {end_x, end_y, ...}) with x-dimension being the inner loop.
        /// </summary>
        /// <template name="integral_t">Integral type for loop iteration</template>
        /// <template name="dimensions">Number of dimensions (= number of loops)</template>
        /// <template name="function_t">Type of function to call with the correct indices</template>
        /// <template name="output_iterator_t">Type of the output iterator which is used to output the return values of the function</template>
        /// <param name="start">First indices</param>
        /// <param name="end">Past the last indices</param>
        /// <param name="function">Function to call with the correct indices</param>
        /// <param name="out_begin">Output iterator pointing at the first element to write the function return value to</param>
        template <typename integral_t, std::size_t dimensions, typename function_t,
            typename output_iterator_t = dummy_output_iterator<typename std::result_of<function_t(Eigen::Matrix<integral_t, dimensions, 1>)>::type>>
        void nested_loop(const Eigen::Matrix<integral_t, dimensions, 1>& start, const Eigen::Matrix<integral_t, dimensions, 1>& end,
            function_t function, output_iterator_t& out_begin = output_iterator_t());

        /// <summary>
        /// Dynamic parallelized nested loop for all permutations in [{start_x, start_y, ...}, {end_x, end_y, ...}) with x-dimension being the inner loop.
        /// </summary>
        /// <template name="integral_t">Integral type for loop iteration</template>
        /// <template name="dimensions">Number of dimensions (= number of loops)</template>
        /// <template name="function_t">Type of function to call with the correct indices</template>
        /// <param name="start">First indices</param>
        /// <param name="end">Past the last indices</param>
        /// <param name="function">Function to call with the correct indices</param>
        /// <param name="out_begin">Output iterator pointing at the first element to write the function return value to (doesn't guarantee sortedness)</param>
        template <typename integral_t, std::size_t dimensions, typename function_t,
            typename output_iterator_t = dummy_output_iterator<typename std::result_of<function_t(Eigen::Matrix<integral_t, dimensions, 1>)>::type>>
        void parallel_nested_loop(const Eigen::Matrix<integral_t, dimensions, 1>& start, const Eigen::Matrix<integral_t, dimensions, 1>& end,
            function_t function, output_iterator_t& out_begin = output_iterator_t());
    }
}

#include "tpf_loop.inl"