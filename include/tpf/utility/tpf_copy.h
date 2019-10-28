#pragma once

#include <memory>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace utility
    {
        /// <summary>
        /// Copy vector of shared pointers using clone
        /// </summary>
        /// <param name="from">Source</param>
        /// <param name="to">Target</param>
        template <typename value_t>
        void copy(const std::vector<std::shared_ptr<value_t>>& from, std::vector<std::shared_ptr<value_t>>& to);

        /// <summary>
        /// Return the provided data unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input data</param>
        /// <returns>Original input data</returns>
        template <typename new_t, typename original_t>
        std::vector<new_t>& copy_data(std::vector<original_t>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

        /// <summary>
        /// Return the provided data unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input data</param>
        /// <returns>Original input data</returns>
        template <typename new_t, typename original_t>
        const std::vector<new_t>& copy_data(const std::vector<original_t>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

        /// <summary>
        /// Copy the data to a new type
        /// </summary>
        /// <template name="new_t">New data type</template>
        /// <template name="original_t">Original data type</template>
        /// <param name="input">Input data</param>
        /// <returns>Copied data</returns>
        template <typename new_t, typename original_t>
        std::vector<new_t> copy_data(const std::vector<original_t>& input, typename std::enable_if<!std::is_same<original_t, new_t>::value>::type* = nullptr);

        /// <summary>
        /// Return the provided grid unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input grid</param>
        /// <returns>Original input grid</returns>
        template <typename new_t, typename original_t>
        std::vector<std::vector<new_t>>& copy_grid(std::vector<std::vector<original_t>>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

        /// <summary>
        /// Return the provided grid unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as NewType</template>
        /// <param name="input">Input grid</param>
        /// <returns>Original input grid</returns>
        template <typename new_t, typename original_t>
        const std::vector<std::vector<new_t>>& copy_grid(const std::vector<std::vector<original_t>>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type* = nullptr);

        /// <summary>
        /// Copy the grid to a new type
        /// </summary>
        /// <template name="new_t">New data type</template>
        /// <template name="original_t">Original data type</template>
        /// <param name="input">Input grid</param>
        /// <returns>Copied grid</returns>
        template <typename new_t, typename original_t>
        std::vector<std::vector<new_t>> copy_grid(const std::vector<std::vector<original_t>>& input, typename std::enable_if<!std::is_same<original_t, new_t>::value>::type* = nullptr);
    }
}

#include "tpf_copy.inl"