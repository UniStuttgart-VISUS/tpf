#include "tpf_copy.h"

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
        inline void copy(const std::vector<std::shared_ptr<value_t>>& from, std::vector<std::shared_ptr<value_t>>& to)
        {
            to.reserve(from.size());

            for (const auto& entry : from)
            {
                to.push_back(entry->clone());
            }
        }

        /// <summary>
        /// Return the provided data unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input data</param>
        /// <returns>Original input data</returns>
        template <typename new_t, typename original_t>
        inline std::vector<new_t>& copy_data(std::vector<original_t>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type*)
        {
            return input;
        }

        /// <summary>
        /// Return the provided data unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input data</param>
        /// <returns>Original input data</returns>
        template <typename new_t, typename original_t>
        inline const std::vector<new_t>& copy_data(const std::vector<original_t>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type*)
        {
            return input;
        }

        /// <summary>
        /// Copy the data to a new type
        /// </summary>
        /// <template name="new_t">New data type</template>
        /// <template name="original_t">Original data type</template>
        /// <param name="input">Input data</param>
        /// <returns>Copied data</returns>
        template <typename new_t, typename original_t>
        inline std::vector<new_t> copy_data(const std::vector<original_t>& input, typename std::enable_if<!std::is_same<original_t, new_t>::value>::type*)
        {
            std::vector<new_t> copy(input.size());

            for (std::size_t i = 0; i < input.size(); ++i)
            {
                copy[i] = static_cast<new_t>(input[i]);
            }

            return copy;
        }

        /// <summary>
        /// Return the provided grid unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as new_t</template>
        /// <param name="input">Input grid</param>
        /// <returns>Original input grid</returns>
        template <typename new_t, typename original_t>
        inline std::vector<std::vector<new_t>>& copy_grid(std::vector<std::vector<original_t>>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type*)
        {
            return input;
        }

        /// <summary>
        /// Return the provided grid unmodified
        /// </summary>
        /// <template name="new_t">Data type</template>
        /// <template name="original_t">Same data type as NewType</template>
        /// <param name="input">Input grid</param>
        /// <returns>Original input grid</returns>
        template <typename new_t, typename original_t>
        inline const std::vector<std::vector<new_t>>& copy_grid(const std::vector<std::vector<original_t>>& input, typename std::enable_if<std::is_same<original_t, new_t>::value>::type*)
        {
            return input;
        }

        /// <summary>
        /// Copy the grid to a new type
        /// </summary>
        /// <template name="new_t">New data type</template>
        /// <template name="original_t">Original data type</template>
        /// <param name="input">Input grid</param>
        /// <returns>Copied grid</returns>
        template <typename new_t, typename original_t>
        inline std::vector<std::vector<new_t>> copy_grid(const std::vector<std::vector<original_t>>& input, typename std::enable_if<!std::is_same<original_t, new_t>::value>::type*)
        {
            std::vector<std::vector<new_t>> copy;

            for (std::size_t d = 0; d < input.size(); ++d)
            {
                copy[d].resize(input[d].size());

                for (std::size_t i = 0; i < input[d].size(); ++i)
                {
                    copy[d][i] = static_cast<new_t>(input[d][i]);
                }
            }

            return copy;
        }
    }
}