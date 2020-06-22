#pragma once

#include "../data/tpf_data_information.h"
#include "../data/tpf_grid.h"
#include "../data/tpf_grid_information.h"

#include "vtkDataArray.h"
#include "vtkRectilinearGrid.h"

#include <string>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        /// <summary>
        /// Load and calculate coordinates, cell sizes and extent
        /// </summary>
        /// <template name="point_t">Point type</template>
        /// <param name="grid">Grid</param>
        /// <param name="data_type">Data type</param>
        /// <param name="cell_coordinates">Output cell coordinates</param>
        /// <param name="node_coordinates">Output node coordinates</param>
        /// <param name="cell_sizes">Output cell sizes</param>
        /// <param name="extent">Output extent</param>
        template <typename point_t>
        void get_grid_information(vtkRectilinearGrid* grid, const data::topology_t data_type, typename data::grid_information<point_t>::array_type& cell_coordinates,
            typename data::grid_information<point_t>::array_type& node_coordinates, typename data::grid_information<point_t>::array_type& cell_sizes, data::extent_t& extent);

        /// <summary>
        /// Save coordinates and extent
        /// </summary>
        /// <template name="point_t">Point type</template>
        /// <param name="grid">Grid</param>
        /// <param name="node_coordinates">Node coordinates</param>
        /// <param name="extent">Extent</param>
        template <typename point_t>
        void set_grid_information(vtkRectilinearGrid* grid, const typename data::grid_information<point_t>::array_type& node_coordinates, const data::extent_t& extent);

        /// <summary>
        /// Create grid from VTK rectilinear grid
        /// </summary>
        /// <template name="data_t">Data type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimensions">Number of dimensions</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        /// <param name="grid">VTK grid</param>
        /// <param name="data_type">Data type</param>
        /// <param name="data_name">Data name</param>
        /// <returns>Grid</returns>
        template <typename data_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns = 1>
        data::grid<data_t, point_t, dimensions, rows, columns> get_grid(vtkRectilinearGrid* grid, data::topology_t data_type, std::string data_name);

        /// <summary>
        /// Create grid from VTK rectilinear grid
        /// </summary>
        /// <template name="data_t">Data type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimensions">Number of dimensions</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        /// <param name="grid">VTK grid</param>
        /// <param name="data_type">Data type</param>
        /// <param name="data_array">Data array extracted from GetInputArrayToProcess</param>
        /// <returns>Grid</returns>
        template <typename data_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns = 1>
        data::grid<data_t, point_t, dimensions, rows, columns> get_grid(vtkRectilinearGrid* grid, data::topology_t data_type, vtkDataArray* data_array);

        /// <summary>
        /// Create grid from VTK rectilinear grid
        /// </summary>
        /// <template name="data_t">Data type</template>
        /// <template name="point_t">Point type</template>
        /// <template name="dimensions">Number of dimensions</template>
        /// <template name="rows">Number of row components</template>
        /// <template name="columns">Number of column components</template>
        /// <param name="grid">VTK grid</param>
        /// <param name="data_type">Data type</param>
        /// <param name="data_names">Data names</param>
        /// <returns>Grid</returns>
        template <typename data_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns = 1>
        data::grid<data_t, point_t, dimensions, rows, columns> get_grid(vtkRectilinearGrid* grid, data::topology_t data_type, const std::vector<std::string>& data_names);
    }
}

#include "tpf_grid.inl"