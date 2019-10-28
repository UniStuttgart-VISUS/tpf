#include "tpf_grid.h"

#include "tpf_data.h"

#include "../data/tpf_data_information.h"
#include "../data/tpf_grid.h"
#include "../data/tpf_grid_information.h"

#include "../log/tpf_log.h"

#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        template <typename point_t>
        inline void get_grid_information(vtkRectilinearGrid* grid, const data::topology_t data_type, typename data::grid_information<point_t>::array_type& cell_coordinates,
            typename data::grid_information<point_t>::array_type& node_coordinates, typename data::grid_information<point_t>::array_type& cell_sizes, data::extent_t& extent)
        {
            try
            {
                // Clear arrays
                extent.clear();
                node_coordinates.clear();
                cell_coordinates.clear();
                cell_sizes.clear();

                extent.reserve(3);
                node_coordinates.reserve(3);
                cell_coordinates.reserve(3);
                cell_sizes.reserve(3);

                // Set extent
                int extent_temp[6];
                grid->GetExtent(extent_temp);

                extent.push_back(std::make_pair(static_cast<std::size_t>(extent_temp[0]), static_cast<std::size_t>(extent_temp[1])));
                extent.push_back(std::make_pair(static_cast<std::size_t>(extent_temp[2]), static_cast<std::size_t>(extent_temp[3])));
                extent.push_back(std::make_pair(static_cast<std::size_t>(extent_temp[4]), static_cast<std::size_t>(extent_temp[5])));

                if (data_type == data::topology_t::CELL_DATA)
                {
                    // Set node coordinates
                    node_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetXCoordinates(), std::string("X coordinates")));
                    node_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetYCoordinates(), std::string("Y coordinates")));
                    node_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetZCoordinates(), std::string("Z coordinates")));

                    // Calculate cell coordinates
                    for (std::size_t c = 0; c < 3; ++c)
                    {
                        cell_coordinates.push_back(std::vector<point_t>(node_coordinates[c].size() - 1));
                        cell_coordinates[c].clear();

                        for (std::size_t i = 0; i < node_coordinates[c].size() - 1; ++i)
                        {
                            cell_coordinates[c].push_back(static_cast<point_t>(0.5L) * (node_coordinates[c][i] + node_coordinates[c][i + 1]));
                        }
                    }

                    --(extent[0].second);
                    --(extent[1].second);
                    --(extent[2].second);
                }
                else
                {
                    // Set cell coordinates
                    cell_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetXCoordinates(), std::string("X coordinates")));
                    cell_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetYCoordinates(), std::string("Y coordinates")));
                    cell_coordinates.push_back(get_data<point_t, vtkFloatArray>(grid->GetZCoordinates(), std::string("Z coordinates")));

                    // Calculate node coordinates
                    for (std::size_t c = 0; c < 3; ++c)
                    {
                        node_coordinates.push_back(std::vector<point_t>(cell_coordinates[c].size() + 1));
                        node_coordinates[c].clear();

                        node_coordinates[c].push_back(static_cast<point_t>(0.0L));

                        for (std::size_t i = 1; i < cell_coordinates[c].size() + 1; ++i)
                        {
                            node_coordinates[c].push_back(static_cast<point_t>(2.0L) * cell_coordinates[c][i - 1] - node_coordinates[c][i - 1]);
                        }
                    }
                }

                // Calculate cell sizes
                for (std::size_t c = 0; c < 3; ++c)
                {
                    cell_sizes.push_back(std::vector<point_t>(cell_coordinates[c].size()));
                    cell_sizes[c].clear();

                    for (std::size_t i = 0; i < cell_coordinates[c].size(); ++i)
                    {
                        cell_sizes[c].push_back(node_coordinates[c][i + 1] - node_coordinates[c][i]);
                    }
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error getting information."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Error getting information."));
            }
        }

        template <typename point_t>
        inline void set_grid_information(vtkRectilinearGrid* grid, const typename data::grid_information<point_t>::array_type& node_coordinates, const data::extent_t& extent)
        {
            // Set extent
            int vtk_extent[6];
            vtk_extent[0] = static_cast<int>(extent[0].first);
            vtk_extent[1] = static_cast<int>(extent[0].second + 1);
            vtk_extent[2] = static_cast<int>(extent[1].first);
            vtk_extent[3] = static_cast<int>(extent[1].second + 1);
            vtk_extent[4] = static_cast<int>(extent[2].first);
            vtk_extent[5] = static_cast<int>(extent[2].second + 1);

            grid->SetExtent(vtk_extent);

            // Set node coordinates
            auto x_coordinates = vtkSmartPointer<vtkDoubleArray>::New();
            auto y_coordinates = vtkSmartPointer<vtkDoubleArray>::New();
            auto z_coordinates = vtkSmartPointer<vtkDoubleArray>::New();

            set_data<point_t, vtkDoubleArray>(x_coordinates.GetPointer(), node_coordinates[0], 1);
            set_data<point_t, vtkDoubleArray>(y_coordinates.GetPointer(), node_coordinates[1], 1);
            set_data<point_t, vtkDoubleArray>(z_coordinates.GetPointer(), node_coordinates[2], 1);

            grid->SetXCoordinates(x_coordinates);
            grid->SetYCoordinates(y_coordinates);
            grid->SetZCoordinates(z_coordinates);
        }

        template <typename data_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline data::grid<data_t, point_t, dimensions, rows, columns> get_grid(vtkRectilinearGrid* grid, const data::topology_t data_type, std::string data_name)
        {
            // Get data
            if (!has_array(grid, data_type, data_name))
            {
                // Try to find an array whose name includes the given name
                if (data_type == data::topology_t::CELL_DATA)
                {
                    for (int i = 0; i < grid->GetCellData()->GetNumberOfArrays(); ++i)
                    {
                        const auto array = grid->GetCellData()->GetAbstractArray(i);

                        if (data_name.find(array->GetName()) != std::string::npos)
                        {
                            data_name = array->GetName();
                            break;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < grid->GetPointData()->GetNumberOfArrays(); ++i)
                    {
                        const auto array = grid->GetPointData()->GetAbstractArray(i);

                        if (data_name.find(array->GetName()) != std::string::npos)
                        {
                            data_name = array->GetName();
                            break;
                        }
                    }
                }

                if (!has_array(grid, data_type, data_name))
                {
                    throw std::runtime_error(__tpf_error_message("Error getting data array. Array with given name not found."));
                }
            }

            auto data = get_data<data_t>(grid, data_type, data_name);

            // Get grid information
            typename data::grid_information<point_t>::array_type cell_coordinates, node_coordinates, cell_sizes;
            data::extent_t extent;

            get_grid_information<point_t>(grid, data::topology_t::CELL_DATA, cell_coordinates, node_coordinates, cell_sizes, extent);

            // Create grid
            data::grid<data_t, point_t, dimensions, rows, columns> output_grid(data_name, extent, std::move(data));
            output_grid.set_grid_information(std::move(cell_coordinates), std::move(node_coordinates), std::move(cell_sizes));

            return output_grid;
        }

        template <typename data_t, typename point_t, std::size_t dimensions, std::size_t rows, std::size_t columns>
        inline data::grid<data_t, point_t, dimensions, rows, columns> get_grid(vtkRectilinearGrid* grid, const data::topology_t data_type, const std::vector<std::string>& data_names)
        {
            // Try to return the grid with the first name
            std::string last_error;

            for (const auto& data_name : data_names)
            {
                try
                {
                    return get_grid<data_t, point_t, dimensions, rows, columns>(grid, data_type, data_name);
                }
                catch (const std::runtime_error& ex)
                {
                    last_error = ex.what();
                }
                catch (...)
                {
                }
            }

            throw std::runtime_error(__tpf_nested_error_message(last_error, "Error getting data array. Array with given names not found."));
        }
    }
}