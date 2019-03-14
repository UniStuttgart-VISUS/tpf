#pragma once

#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// Topology type
        enum class topology_t {
            CELL_DATA, POINT_DATA, OBJECT_DATA
        };

        /// <summary>
        /// Struct specifying a data set
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        template <typename value_t, std::size_t rows = 1, std::size_t columns = 1>
        struct data_information
        {
            /// Template information
            using value_type = value_t;

            constexpr std::size_t get_num_rows() { return rows; }
            constexpr std::size_t get_num_columns() { return columns; }

            /// Data set name and topology
            std::string name;
            topology_t topology;
        };

        /// Convenience for (spatial) data information
        using extent_t = std::vector<std::pair<std::size_t, std::size_t>>;
        template <typename arithmetic_t> using bbox_t = std::vector<std::pair<arithmetic_t, arithmetic_t>>;
        template <typename arithmetic_t> using area_t = std::vector<std::pair<arithmetic_t, arithmetic_t>>;
    }
}