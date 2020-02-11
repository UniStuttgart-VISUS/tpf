#include "tpf_polydata_element.h"

#include "Eigen/Dense"

#include <string>
#include <vector>

namespace tpf
{
    namespace data
    {
        template<typename value_t, std::size_t rows, std::size_t columns>
        inline polydata_element<value_t, rows, columns>::polydata_element(const std::string& name, topology_t topology, const std::vector<Eigen::Matrix<value_t, rows, columns>>& values)
            : name(name), topology(topology), values(values)
        { }

        template<typename value_t>
        inline polydata_element<value_t, 1, 1>::polydata_element(const std::string& name, topology_t topology, const std::vector<value_t>& values)
            : name(name), topology(topology), values(values)
        { }
    }
}