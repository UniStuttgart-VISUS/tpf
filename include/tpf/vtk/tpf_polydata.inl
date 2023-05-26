#include "tpf_polydata.h"

#include "tpf_data.h"

#include "../data/tpf_array.h"
#include "../data/tpf_data_information.h"
#include "../data/tpf_polydata.h"

#include "../geometry/tpf_geometric_object.h"
#include "../geometry/tpf_line.h"
#include "../geometry/tpf_point.h"
#include "../geometry/tpf_polygon.h"

#include "../log/tpf_log.h"

#include "../math/tpf_vector.h"

#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <algorithm>
#include <exception>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        namespace 
        {
            /// <summary>
            /// Mesh storing points, vertices, lines and poligonal cells
            /// </summary>
            struct mesh
            {
                vtkSmartPointer<vtkPoints> points;
                vtkSmartPointer<vtkCellArray> vertices;
                vtkSmartPointer<vtkCellArray> lines;
                vtkSmartPointer<vtkCellArray> polygonals;
            };

            /// <summary>
            /// Create polygonal mesh from geometric objects
            /// </summary>
            /// <param name="objects">Geometric objects</param>
            /// <returns>Points and cells</returns>
            template <typename float_t>
            inline mesh create_mesh(const std::vector<std::shared_ptr<geometry::geometric_object<float_t>>>& objects)
            {
                try
                {
                    const auto is_vertex = static_cast<uint64_t>(geometry::geometry_t::POINT);
                    const auto is_line =
                        static_cast<uint64_t>(geometry::geometry_t::LINE) |
                        static_cast<uint64_t>(geometry::geometry_t::LINE_STRIP);
                    const auto is_polygon =
                        static_cast<uint64_t>(geometry::geometry_t::TRIANGLE) |
                        static_cast<uint64_t>(geometry::geometry_t::TETRAHEDRON) |
                        static_cast<uint64_t>(geometry::geometry_t::RECTANGLE) |
                        static_cast<uint64_t>(geometry::geometry_t::CUBOID) |
                        static_cast<uint64_t>(geometry::geometry_t::POLYGON) |
                        static_cast<uint64_t>(geometry::geometry_t::POLYHEDRON) |
                        static_cast<uint64_t>(geometry::geometry_t::MESH);

                    // Get number of points and cell indices
                    std::size_t num_points = 0;
                    std::size_t num_verts = 0;
                    std::size_t num_vert_indices = 0;
                    std::size_t num_lines = 0;
                    std::size_t num_line_indices = 0;
                    std::size_t num_polys = 0;
                    std::size_t num_polys_indices = 0;

                    std::for_each(objects.begin(), objects.end(), [&num_points](const std::shared_ptr<geometry::geometric_object<float_t>>& obj) { num_points += obj->get_num_points(); });
                    std::for_each(objects.begin(), objects.end(), [&is_vertex, &is_line, &is_polygon, &num_verts, &num_vert_indices, &num_lines,
                        &num_line_indices, &num_polys, &num_polys_indices](const std::shared_ptr<geometry::geometric_object<float_t>>& obj)
                    {
                        if (static_cast<uint64_t>(obj->get_type()) & is_vertex)
                        {
                            ++num_verts;
                            num_vert_indices += 2;
                        }
                        else if (static_cast<uint64_t>(obj->get_type()) & is_line)
                        {
                            ++num_lines;
                            num_line_indices += obj->get_cells()[0].size() + 1;
                        }
                        else if (static_cast<uint64_t>(obj->get_type()) & is_polygon)
                        {
                            for (const auto& cell : obj->get_cells())
                            {
                                ++num_polys;
                                num_polys_indices += cell.size() + 1;
                            }
                        }
                    });

                    // Create points and cell indices
                    auto points = vtkSmartPointer<vtkPoints>::New();
                    points->SetNumberOfPoints(static_cast<int>(num_points));

                    auto vert_indices = vtkSmartPointer<vtkIdTypeArray>::New();
                    vert_indices->SetNumberOfComponents(1);
                    vert_indices->SetNumberOfTuples(static_cast<int>(num_vert_indices));

                    auto line_indices = vtkSmartPointer<vtkIdTypeArray>::New();
                    line_indices->SetNumberOfComponents(1);
                    line_indices->SetNumberOfTuples(static_cast<int>(num_line_indices));

                    auto poly_indices = vtkSmartPointer<vtkIdTypeArray>::New();
                    poly_indices->SetNumberOfComponents(1);
                    poly_indices->SetNumberOfTuples(static_cast<int>(num_polys_indices));

                    std::size_t point_index = 0;
                    std::size_t vert_index = 0;
                    std::size_t line_index = 0;
                    std::size_t poly_index = 0;

                    for (const std::shared_ptr<geometry::geometric_object<float_t>> obj : objects)
                    {
                        const auto obj_points = obj->get_points();

                        for (const auto& cell : obj->get_cells())
                        {
                            if (static_cast<uint64_t>(obj->get_type()) & is_vertex)
                            {
                                vert_indices->SetValue(static_cast<int>(vert_index++), static_cast<vtkIdType>(1));
                            }
                            else if (static_cast<uint64_t>(obj->get_type()) & is_line)
                            {
                                line_indices->SetValue(static_cast<int>(line_index++), static_cast<vtkIdType>(cell.size()));
                            }
                            else if (static_cast<uint64_t>(obj->get_type()) & is_polygon)
                            {
                                poly_indices->SetValue(static_cast<int>(poly_index++), static_cast<vtkIdType>(cell.size()));
                            }

                            for (const auto& index : cell)
                            {
                                float p[3] = { static_cast<float>(obj_points[index][0]), static_cast<float>(obj_points[index][1]), static_cast<float>(obj_points[index][2]) };
                                points->SetPoint(static_cast<vtkIdType>(point_index + index), p);

                                if (static_cast<uint64_t>(obj->get_type()) & is_vertex)
                                {
                                    vert_indices->SetValue(static_cast<int>(vert_index++), static_cast<vtkIdType>(point_index + index));
                                }
                                else if (static_cast<uint64_t>(obj->get_type()) & is_line)
                                {
                                    line_indices->SetValue(static_cast<int>(line_index++), static_cast<vtkIdType>(point_index + index));
                                }
                                else if (static_cast<uint64_t>(obj->get_type()) & is_polygon)
                                {
                                    poly_indices->SetValue(static_cast<int>(poly_index++), static_cast<vtkIdType>(point_index + index));
                                }
                            }
                        }

                        point_index += obj_points.size();
                    }

                    // Create vertices
                    vtkSmartPointer<vtkCellArray> verts = nullptr;

                    if (num_verts > 0)
                    {
                        verts = vtkSmartPointer<vtkCellArray>::New();
                        verts->SetCells(static_cast<vtkIdType>(num_verts), vert_indices);
                    }

                    // Create lines
                    vtkSmartPointer<vtkCellArray> lines = nullptr;

                    if (num_lines > 0)
                    {
                        lines = vtkSmartPointer<vtkCellArray>::New();
                        lines->SetCells(static_cast<vtkIdType>(num_lines), line_indices);
                    }

                    // Create polygons
                    vtkSmartPointer<vtkCellArray> polys = nullptr;

                    if (num_polys > 0)
                    {
                        polys = vtkSmartPointer<vtkCellArray>::New();
                        polys->SetCells(static_cast<vtkIdType>(num_polys), poly_indices);
                    }

                    return{ points, verts, lines, polys };
                }
                catch (const std::exception& ex)
                {
                    throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error creating mesh."));
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error creating mesh."));
                }
            }

            /// <summary>
            /// Create geometric objects from polygonal mesh
            /// </summary>
            /// <param name="input_mesh">Input mesh</param>
            /// <returns>Geometric data</returns>
            template <typename float_t>
            inline std::vector<std::shared_ptr<geometry::geometric_object<float_t>>> create_geometry(mesh input_mesh)
            {
                try
                {
                    std::vector<std::shared_ptr<geometry::geometric_object<float_t>>> geometric_objects;

                    auto points = input_mesh.points;
                    auto verts = input_mesh.vertices;
                    auto lines = input_mesh.lines;
                    auto polys = input_mesh.polygonals;

                    auto vtk_point_ids = vtkSmartPointer<vtkIdList>::New();
                    double vtk_point[3];

                    // Create points (vertices)
                    verts->InitTraversal();

                    while (verts->GetNextCell(vtk_point_ids))
                    {
                        points->GetPoint(vtk_point_ids->GetId(0), vtk_point);
                        const auto point_1 = math::vec3_t<float_t>(static_cast<float_t>(vtk_point[0]), static_cast<float_t>(vtk_point[1]), static_cast<float_t>(vtk_point[2]));

                        geometric_objects.push_back(std::make_shared<geometry::point<float_t>>(point_1));
                    }

                    // Create lines
                    lines->InitTraversal();

                    while (lines->GetNextCell(vtk_point_ids))
                    {
                        points->GetPoint(vtk_point_ids->GetId(0), vtk_point);
                        const auto point_1 = math::vec3_t<float_t>(static_cast<float_t>(vtk_point[0]), static_cast<float_t>(vtk_point[1]), static_cast<float_t>(vtk_point[2]));

                        points->GetPoint(vtk_point_ids->GetId(1), vtk_point);
                        const auto point_2 = math::vec3_t<float_t>(static_cast<float_t>(vtk_point[0]), static_cast<float_t>(vtk_point[1]), static_cast<float_t>(vtk_point[2]));

                        geometric_objects.push_back(std::make_shared<geometry::line<float_t>>(point_1, point_2));
                    }

                    // Create polygonals
                    polys->InitTraversal();

                    while (polys->GetNextCell(vtk_point_ids))
                    {
                        std::vector<geometry::point<float_t>> polygon_vertices;
                        polygon_vertices.reserve(vtk_point_ids->GetNumberOfIds());

                        for (vtkIdType i = 0; i < vtk_point_ids->GetNumberOfIds(); ++i)
                        {
                            points->GetPoint(vtk_point_ids->GetId(i), vtk_point);

                            polygon_vertices.emplace_back(static_cast<float_t>(vtk_point[0]), static_cast<float_t>(vtk_point[1]), static_cast<float_t>(vtk_point[2]));
                        }

                        geometric_objects.push_back(std::make_shared<geometry::polygon<float_t>>(polygon_vertices, true, false));
                    }

                    return geometric_objects;
                }
                catch (const std::exception& ex)
                {
                    throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Error creating geometry."));
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error creating geometry."));
                }
            }

            template <typename point_t>
            inline void get_polydata_array(vtkPolyData* data, data::polydata<point_t>& polydata)
            { }

            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns>
            inline void get_polydata_array(vtkPolyData* data, data::polydata<point_t>& polydata,
                const data::data_information<value_t, rows, columns>& data_info)
            {
                if (!has_array(data, data_info.topology, data_info.name))
                {
                    throw std::runtime_error(__tpf_error_message("Error getting data array. Array with given name not found."));
                }

                polydata.add(std::make_shared<data::array<value_t, rows, columns>>(data_info.name, std::move(get_data<value_t>(data, data_info.topology, data_info.name))), data_info.topology);
            }

            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_info_t>
            inline void get_polydata_array(vtkPolyData* data, data::polydata<point_t>& polydata,
                const data::data_information<value_t, rows, columns>& data_info, data_info_t... rest)
            {
                get_polydata_array(data, polydata, data_info);
                get_polydata_array(data, polydata, rest...);
            }

            template <typename point_t>
            inline void set_polydata_array(vtkPolyData* data, const data::polydata<point_t>& polydata)
            { }

            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_info_t>
            inline void set_polydata_array(vtkPolyData* data, const data::polydata<point_t>& polydata,
                const data::data_information<value_t, rows, columns>& data_info)
            {
                data::topology_t topology = data_info.topology;
                std::shared_ptr<data::array<value_t, rows, columns>> data_array = nullptr;

                if (data_info.topology == data::topology_t::CELL_DATA)
                {
                    if (!polydata.has_cell_data(data_info.name))
                    {
                        throw std::runtime_error(__tpf_error_message("Error getting cell data array. Array with given name not found."));
                    }

                    data_array = polydata.template get_cell_data_as<value_t, rows, columns>(data_info.name);
                }
                else if (data_info.topology == data::topology_t::POINT_DATA || data_info.topology == data::topology_t::TEXTURE_COORDINATES)
                {
                    if (!polydata.has_point_data(data_info.name))
                    {
                        throw std::runtime_error(__tpf_error_message("Error getting point data array. Array with given name not found."));
                    }

                    data_array = polydata.template get_point_data_as<value_t, rows, columns>(data_info.name);
                }
                else if (data_info.topology == data::topology_t::OBJECT_DATA)
                {
                    if (!polydata.has_object_data(data_info.name))
                    {
                        throw std::runtime_error(__tpf_error_message("Error getting object data array. Array with given name not found."));
                    }

                    // Map object data to cell data.
                    const auto object_data_array = polydata.template get_object_data_as<value_t, rows, columns>(data_info.name);
                    const auto cell_nums = polydata.get_num_object_cells();

                    if (object_data_array != nullptr)
                    {
                        topology = data::topology_t::CELL_DATA;
                        data_array = std::make_shared<data::array<value_t, rows, columns>>(data_info.name, polydata.get_num_cells());

                        std::size_t idx = 0;
                        for (std::size_t i = 0; i < polydata.get_num_objects(); ++i)
                        {
                            for (std::size_t j = 0; j < cell_nums[i]; ++j)
                            {
                                data_array->at(idx) = object_data_array->at(i);
                                ++idx;
                            }
                        }
                    }
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Invalid topology type."));
                }

                if (data_array == nullptr)
                {
                    throw std::runtime_error(__tpf_error_message("Data array found but incompatible with the given data type, and number of rows and columns."));
                }

                set_data<value_t>(data, topology, data_info.name, data_array->get_data(), data_array->get_num_components());
            }

            template <typename point_t, typename value_t, std::size_t rows, std::size_t columns, typename... data_info_t>
            inline void set_polydata_array(vtkPolyData* data, const data::polydata<point_t>& polydata,
                const data::data_information<value_t, rows, columns>& data_info, const data_info_t&... rest)
            {
                set_polydata_array(data, polydata, data_info);
                set_polydata_array(data, polydata, rest...);
            }
        }

        template <typename point_t, typename... data_info_t>
        inline data::polydata<point_t> get_polydata(vtkPolyData* data, const data_info_t&... data_info)
        {
            if (data == nullptr)
            {
                return data::polydata<point_t>();
            }

            mesh input_mesh = { vtkSmartPointer<vtkPoints>(data->GetPoints()), vtkSmartPointer<vtkCellArray>(data->GetVerts()),
                vtkSmartPointer<vtkCellArray>(data->GetLines()), vtkSmartPointer<vtkCellArray>(data->GetPolys()) };

            auto geometry = create_geometry<point_t>(input_mesh);
            auto polydata = data::polydata<point_t>(geometry);

            get_polydata_array<point_t>(data, polydata, data_info...);

            return polydata;
        }

        template <typename point_t, typename... data_info_t>
        inline void set_polydata(vtkPolyData* data, const data::polydata<point_t>& polydata, const data_info_t&... data_info)
        {
            mesh mesh = create_mesh<point_t>(polydata.get_geometry());

            set_polydata_array<point_t>(data, polydata, data_info...);

            data->SetPoints(mesh.points);
            data->SetVerts(mesh.vertices);
            data->SetLines(mesh.lines);
            data->SetPolys(mesh.polygonals);
        }

        template <typename point_t, typename... data_info_t>
        inline void append_polydata(vtkPolyData* data, const data::polydata<point_t>& polydata, const data_info_t&... data_info)
        {
            set_polydata_array<point_t>(data, polydata, data_info...);
        }
    }
}
