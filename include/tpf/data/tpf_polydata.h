#pragma once

#include "tpf_array.h"
#include "tpf_array_base.h"
#include "tpf_data_information.h"
#include "tpf_polydata_element.h"

#include "../geometry/tpf_geometric_object.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Class for storing polydata and additional data arrays for point or cell data
        /// </summary>
        /// <template name="point_t">Data type for storing positional data</template>
        template <typename point_t>
        class polydata
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            polydata();

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="geometry">Geometry</param>
            polydata(const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& geometry);

            /// <summary>
            /// Copy constructor
            /// </summary>
            /// <param name="copy">Copy</param>
            polydata(const polydata& copy);

            /// <summary>
            /// Move constructor
            /// </summary>
            /// <param name="move">Move</param>
            polydata(polydata&& move);

            /// <summary>
            /// Copy assignment
            /// </summary>
            /// <param name="copy">Copy</param>
            /// <returns>This</returns>
            polydata& operator=(const polydata& copy);

            /// <summary>
            /// Move assignment
            /// </summary>
            /// <param name="move">Move</param>
            /// <returns>This</returns>
            polydata& operator=(polydata&& move);

            /// <summary>
            /// Insert geometric object without additional data (only allowed if no data array has yet been appended)
            /// </summary>
            /// <param name="object">Geometric object</param>
            /// <returns>Index of inserted object</returns>
            std::size_t insert(std::shared_ptr<geometry::geometric_object<point_t>> object);

            /// <summary>
            /// Insert geometric objects without additional data (only allowed if no data array has yet been appended)
            /// </summary>
            /// <param name="objects">Geometric objects</param>
            /// <returns>Index of last inserted object</returns>
            std::size_t insert(std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& objects);

            /// <summary>
            /// Insert geometric object with additional data
            /// </summary>
            /// <template name="value_t">Value type of first attached data element</template>
            /// <template name="rows">Number of rows of first attached data element</template>
            /// <template name="columns">Number of columns of first attached data element</template>
            /// <template name="data_elements_t">Types of further attached data elements</template>
            /// <param name="object">Geometric object</param>
            /// <param name="first_element">First attached data element</param>
            /// <param name="data_elements">Further attached data elements</param>
            /// <returns>Index of inserted object</returns>
            template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
            std::size_t insert(std::shared_ptr<geometry::geometric_object<point_t>> object,
                const polydata_element<value_t, rows, columns>& first_element, const data_elements_t&... data_elements);

            /// <summary>
            /// Add point or cell data array
            /// </summary>
            /// <param name="array">Data array</param>
            /// <param name="type">Point or cell data</param>
            void add(std::shared_ptr<array_base> array, topology_t type);

            /// <summary>
            /// Merge another poly data object into this one
            /// </summary>
            /// <param name="other">Other poly data</param>
            void merge(const polydata& other);

            /// <summary>
            /// Create and add point or cell data array
            /// </summary>
            /// <template name="value_t">Value type</template>
            /// <template name="rows">Number of row components</template>
            /// <template name="columns">Number of column components</template>
            /// <param name="name">Data array name</param>
            /// <param name="type">Point or cell data</param>
            /// <returns>New, added data array</returns>
            template <typename value_t, std::size_t rows, std::size_t columns = 1>
            std::shared_ptr<array<value_t, rows, columns>> create(const std::string& name, topology_t type);

            /// <summary>
            /// Extract an area, given by the extent
            /// </summary>
            /// <template name="AltPT">Alternative point type for the extracted area</template>
            /// <param name="subarea">Area to extract</param>
            /// <param name="include_arrays">Decide to also extract data arrays</param>
            /// <returns>Polydata</returns>
            polydata<point_t> extract_area(const area_t<point_t>& subarea, bool include_arrays = false) const;

            polydata<point_t> extract_area(const area_t<point_t>& subarea, bool include_arrays = false);

            /// <summary>
            /// Extract specific arrays of an area, given by the extent
            /// </summary>
            /// <template name="AltPT">Alternative point type for the extracted area</template>
            /// <template name="data_elements_t">Types of attached data elements</template>
            /// <param name="subarea">Area to extract</param>
            /// <param name="data_elements">Attached data elements</param>
            /// <returns>Polydata</returns>
            template <typename... data_elements_t>
            polydata<point_t> extract_area(const area_t<point_t>& subarea, const data_elements_t&... data_elements) const;

            template <typename... data_elements_t>
            polydata<point_t> extract_area(const area_t<point_t>& subarea, const data_elements_t&... data_elements);

            /// <summary>
            /// Get geometric object by its ID
            /// </summary>
            /// <param name="id">ID</param>
            /// <returns>Geometric object</returns>
            const std::shared_ptr<geometry::geometric_object<point_t>> get_object(std::size_t id) const;

            /// <summary>
            /// Return geometric objects
            /// </summary>
            /// <returns>Geometric objects</returns>
            const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& get_geometry() const;

            /// <summary>
            /// Check for the existance of a certain point data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>True if it exists</returns>
            bool has_point_data(const std::string& name) const;

            /// <summary>
            /// Return point data arrays
            /// </summary>
            /// <returns>Cell point arrays</returns>
            const std::vector<std::shared_ptr<array_base>>& get_point_data() const;

            /// <summary>
            /// Return requested point data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>Requested array, nullptr if not found</returns>
            const std::shared_ptr<array_base> get_point_data(const std::string& name) const;

            /// <summary>
            /// Return requested point data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <template name="value_t">Value type</template>
            /// <template name="rows">Number of row components</template>
            /// <template name="columns">Number of column components</template>
            /// <returns>Requested array, nullptr if not found</returns>
            template <typename value_t, std::size_t rows, std::size_t columns = 1>
            const std::shared_ptr<array<value_t, rows, columns>> get_point_data_as(const std::string& name) const;

            /// <summary>
            /// Check for the existance of a certain cell data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>True if it exists</returns>
            bool has_cell_data(const std::string& name) const;

            /// <summary>
            /// Return cell data arrays
            /// </summary>
            /// <returns>Cell data arrays</returns>
            const std::vector<std::shared_ptr<array_base>>& get_cell_data() const;

            /// <summary>
            /// Return requested cell data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>Requested array, nullptr if not found</returns>
            const std::shared_ptr<array_base> get_cell_data(const std::string& name) const;

            /// <summary>
            /// Return requested cell data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <template name="value_t">Value type</template>
            /// <template name="rows">Number of row components</template>
            /// <template name="columns">Number of column components</template>
            /// <returns>Requested array, nullptr if not found</returns>
            template <typename value_t, std::size_t rows, std::size_t columns = 1>
            const std::shared_ptr<array<value_t, rows, columns>> get_cell_data_as(const std::string& name) const;

            /// <summary>
            /// Check for the existance of a certain object data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>True if it exists</returns>
            bool has_object_data(const std::string& name) const;

            /// <summary>
            /// Return object data arrays
            /// </summary>
            /// <returns>Object data arrays</returns>
            const std::vector<std::shared_ptr<array_base>>& get_object_data() const;

            /// <summary>
            /// Return requested object data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <returns>Requested array, nullptr if not found</returns>
            const std::shared_ptr<array_base> get_object_data(const std::string& name) const;

            /// <summary>
            /// Return requested object data array
            /// </summary>
            /// <param name="name">Array name</param>
            /// <template name="value_t">Value type</template>
            /// <template name="rows">Number of row components</template>
            /// <template name="columns">Number of column components</template>
            /// <returns>Requested array, nullptr if not found</returns>
            template <typename value_t, std::size_t rows, std::size_t columns = 1>
            const std::shared_ptr<array<value_t, rows, columns>> get_object_data_as(const std::string& name) const;

            /// <summary>
            /// Return number of points
            /// </summary>
            /// <returns>Number of points</returns>
            std::size_t get_num_points() const;

            /// <summary>
            /// Return number of cells
            /// </summary>
            /// <returns>Number of cells</returns>
            std::size_t get_num_cells() const;

            /// <summary>
            /// Return number of objects
            /// </summary>
            /// <returns>Number of objects</returns>
            std::size_t get_num_objects() const;

            /// <summary>
            /// Return list of number of cells per object
            /// </summary>
            /// <returns>Number of cells per objects</returns>
            std::vector<std::size_t> get_num_object_cells() const;

        private:
            /// <summary>
            /// Insert geometric object with additional data
            /// </summary>
            /// <template name="value_t">Value type of first attached data element</template>
            /// <template name="rows">Number of rows of first attached data element</template>
            /// <template name="columns">Number of columns of first attached data element</template>
            /// <template name="data_elements_t">Types of further attached data elements</template>
            /// <param name="first_element">First attached data element</param>
            /// <param name="data_elements">Further attached data elements</param>
            /// <returns>Number of attached data elements</returns>
            template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
            std::size_t insert(std::size_t num_points, std::size_t num_cells, const polydata_element<value_t, rows, columns>& first_element, const data_elements_t&... data_elements);

            /// <summary>
            /// Insert geometric object with additional data
            /// </summary>
            /// <template name="value_t">Value type of first attached data element</template>
            /// <template name="rows">Number of rows of first attached data element</template>
            /// <template name="columns">Number of columns of first attached data element</template>
            /// <param name="element">Attached data element</param>
            /// <returns>Number of attached data elements</returns>
            template <typename value_t, std::size_t rows, std::size_t columns>
            std::size_t insert(std::size_t num_points, std::size_t num_cells, const polydata_element<value_t, rows, columns>& element);

            /// <summary>
            /// Extract indices corresponding to the sub extent given
            /// </summary>
            /// <param name="subarea">Area to extract</param>
            /// <returns>Indices of objects and array elements to be extracted</returns>
            std::vector<std::size_t> extract_indices(const area_t<point_t>& subarea) const;

            /// <summary>
            /// Extract the geometries, given by the indices
            /// </summary>
            /// <param name="indices">Indices of the geometry to extract</param>
            /// <returns>Geometry</returns>
            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> extract_geometry(const std::vector<std::size_t>& indices) const;

            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> extract_geometry(const std::vector<std::size_t>& indices);

            /// <summary>
            ///  Extract specific arrays, with elements given by the indices
            /// </summary>
            /// <param name="indices">Indices of the elements to extract</param>
            /// <param name="first_element">First data array</param>
            /// <param name="data_elements">Further data arrays</param>
            /// <template name="value_t">Value type of the first array</template>
            /// <template name="rows">Number of rows of elements in the first array</template>
            /// <template name="columns">Number of columns of elements in the first array</template>
            /// <template name="data_elements_t">Types of further data information objects</template>
            /// <returns>Arrays extracted</returns>
            template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
            std::vector<std::pair<std::shared_ptr<data::array_base>, data::topology_t>> extract_array(const std::vector<std::size_t>& indices,
                const data_information<value_t, rows, columns>& first_element, const data_elements_t&... data_elements) const;

            /// <summary>
            ///  Extract specific array, with elements given by the indices
            /// </summary>
            /// <param name="indices">Indices of the elements to extract</param>
            /// <param name="element">Data array</param>
            /// <template name="value_t">Value type of the array</template>
            /// <template name="rows">Number of rows of elements in the array</template>
            /// <template name="columns">Number of columns of elements in the array</template>
            /// <returns>Array extracted</returns>
            template <typename value_t, std::size_t rows, std::size_t columns>
            std::vector<std::pair<std::shared_ptr<data::array_base>, data::topology_t>> extract_array(
                const std::vector<std::size_t>& indices, const data_information<value_t, rows, columns>& element) const;

            /// Geometric data
            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> geometry;

            /// Additional data arrays
            std::vector<std::shared_ptr<array_base>> point_data;
            std::vector<std::shared_ptr<array_base>> cell_data;
            std::vector<std::shared_ptr<array_base>> object_data;
        };
    }
}

#include "tpf_polydata.inl"
