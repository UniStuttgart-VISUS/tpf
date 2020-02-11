#pragma once

#include "Eigen/Dense"

#include <string>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Struct for a vector data element
        /// </summary>
        /// <template name="value_t">Value type</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        template <typename value_t, std::size_t rows, std::size_t columns = 1>
        struct polydata_element
        {
            /// Template parameters
            using value_type = Eigen::Matrix<value_t, rows, columns>;
            const std::size_t num_rows = rows;
            const std::size_t num_columns = columns;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="topology">Topology type</param>
            /// <param name="values">Value(s)</param>
            polydata_element(const std::string& name, topology_t topology, const std::vector<Eigen::Matrix<value_t, rows, columns>>& values);

            /// Data set name and topology
            const std::string name;
            const topology_t topology;

            /// Data values
            const std::vector<Eigen::Matrix<value_t, rows, columns>>& values;
        };

        /// <summary>
        /// Struct for a scalar data element
        /// </summary>
        /// <template name="value_t">Value type</template>
        template <typename value_t>
        struct polydata_element<value_t, 1, 1>
        {
            /// Template parameters
            using value_type = value_t;
            const std::size_t num_rows = 1;
            const std::size_t num_columns = 1;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="topology">Topology type</param>
            /// <param name="values">Value(s)</param>
            polydata_element(const std::string& name, topology_t topology, const std::vector<value_t>& values);

            /// Data set name and topology
            const std::string name;
            const topology_t topology;

            /// Data values
            const std::vector<value_t>& values;
        };
    }
}

#include "tpf_polydata_element.inl"