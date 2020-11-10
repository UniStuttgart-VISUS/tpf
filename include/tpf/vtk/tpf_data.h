#pragma once

#include "tpf_vtk_traits.h"

#include "../data/tpf_array.h"
#include "../data/tpf_data_information.h"

#include "vtkAbstractArray.h"
#include "vtkDataSet.h"

#include <string>
#include <vector>

#ifdef _MSC_EXTENSIONS
#include <iterator>
#include <type_traits>
#endif

namespace tpf
{
    namespace vtk
    {
#ifdef _MSC_EXTENSIONS
#define CHECKED_ITERATOR(x, array_t) stdext::checked_array_iterator<decltype(x->Begin())>(x->Begin(), x->GetNumberOfValues())
#else
#define CHECKED_ITERATOR(x, array_t) x->Begin()
#endif

        /// <summary>
        /// Warn about a mismatch in types
        /// </summary>
        /// <param name="array_name">Array name</param>
        template <typename expected_t, typename found_t>
        void warn_mismatch(const std::string& array_name);

        /// <summary>
        /// Does the data set contain any data?
        /// </summary>
        /// <param name="data">Data set</param>
        /// <param name="data_type">Data type</param>
        /// <return>Returns true if data set contains data, false otherwise</return>
        bool has_data(vtkDataSet* data, const data::topology_t data_type);

        /// <summary>
        /// Does the data set contain the specified array?
        /// </summary>
        /// <param name="data">Data set</param>
        /// <param name="data_type">Data type</param>
        /// <param name="array_name">Array name</param>
        /// <return>Returns true if data set contains the specified array, false otherwise</return>
        bool has_array(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name);

        /// <summary>
        /// Load data from data array
        /// </summary>
        /// <template name="value_t">Scalar type</template>
        /// <template name="array_t">Derived class of vtkDataArray</template>
        /// <template name="dynamic">Allow dynamic downcast if array type doesn't match</template>
        /// <param name="data">Data array</param>
        /// <param name="array_name">Array name</param>
        /// <return>Array</return>
        template <typename value_t, class array_t = typename vtk_array<value_t>::type, bool dynamic = true>
        std::vector<value_t> get_data(vtkAbstractArray* data, const std::string& array_name);

        /// <summary>
        /// Load data from data set
        /// </summary>
        /// <template name="value_t">Scalar type</template>
        /// <template name="array_t">Derived class of vtkDataArray</template>
        /// <template name="dynamic">Allow dynamic downcast if array type doesn't match</template>
        /// <param name="data">Data set</param>
        /// <param name="data_type">Data type</param>
        /// <param name="array_name">Array name</param>
        /// <return>Array</return>
        template <typename value_t, class array_t = typename vtk_array<value_t>::type, bool dynamic = true>
        std::vector<value_t> get_data(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name);

        /// <summary>
        /// Save data array to data set
        /// </summary>
        /// <template name="value_t">Scalar type</template>
        /// <template name="array_t">Derived class of vtkDataArray</template>
        /// <param name="data">Data array</param>
        /// <param name="data_array">Data array to save</param>
        /// <param name="num_components">Number of components in the data array</param>
        template <typename value_t, class array_t = typename vtk_array<value_t>::type>
        void set_data(vtkAbstractArray* data, const std::vector<value_t>& data_array, const std::size_t num_components);

        /// <summary>
        /// Save data array to data set
        /// </summary>
        /// <template name="value_t">Scalar type</template>
        /// <template name="array_t">Derived class of vtkDataArray</template>
        /// <param name="data">Data set</param>
        /// <param name="data_type">Data type</param>
        /// <param name="array_name">Array name</param>
        /// <param name="data_array">Data array to save</param>
        /// <param name="num_components">Number of components in the data array</param>
        /// <return>Array</return>
        template <typename value_t, class array_t = typename vtk_array<value_t>::type>
        void set_data(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name,
            const std::vector<value_t>& data_array, const std::size_t num_components);

        /// <summary>
        /// Save data array to data set
        /// </summary>
        /// <template name="value_t">Scalar type</template>
        /// <template name="array_t">Derived class of vtkDataArray</template>
        /// <param name="data">Data set</param>
        /// <param name="data_type">Data type</param>
        /// <param name="array">Array to save</param>
        /// <return>Array</return>
        template <typename value_t, class array_t = typename vtk_array<value_t>::type, int rows = 1, int columns = 1>
        void set_data(vtkDataSet* data, const data::topology_t data_type, const data::array<value_t, rows, columns>& array);
    }
}

#include "tpf_data.inl"
