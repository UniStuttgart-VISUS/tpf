#include "tpf_polydata.h"

#include "tpf_array.h"
#include "tpf_array_base.h"
#include "tpf_data_information.h"
#include "tpf_polydata_element.h"

#include "../geometry/tpf_geometric_object.h"

#include "../log/tpf_log.h"

#include "../utility/tpf_copy.h"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename point_t>
        inline polydata<point_t>::polydata() { }

        template <typename point_t>
        inline polydata<point_t>::polydata(const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& geometry)
        {
            this->geometry = geometry;
        }

        template <typename point_t>
        inline polydata<point_t>::polydata(const polydata& copy)
        {
            utility::copy(copy.geometry, this->geometry);

            utility::copy(copy.point_data, this->point_data);
            utility::copy(copy.cell_data, this->cell_data);
            utility::copy(copy.object_data, this->object_data);
        }

        template <typename point_t>
        inline polydata<point_t>::polydata(polydata&& move)
        {
            this->geometry = move.geometry;
            move.geometry.clear();

            this->point_data = move.point_data;
            move.point_data.clear();

            this->cell_data = move.cell_data;
            move.cell_data.clear();

            this->object_data = move.object_data;
            move.object_data.clear();
        }

        template <typename point_t>
        inline polydata<point_t>& polydata<point_t>::operator=(const polydata& copy)
        {
            utility::copy(copy.geometry, this->geometry);

            utility::copy(copy.point_data, this->point_data);
            utility::copy(copy.cell_data, this->cell_data);
            utility::copy(copy.object_data, this->object_data);

            return *this;
        }

        template <typename point_t>
        inline polydata<point_t>& polydata<point_t>::operator=(polydata&& move)
        {
            this->geometry = move.geometry;
            move.geometry.clear();

            this->point_data = move.point_data;
            move.point_data.clear();

            this->cell_data = move.cell_data;
            move.cell_data.clear();

            this->object_data = move.object_data;
            move.object_data.clear();

            return *this;
        }

        template <typename point_t>
        inline std::size_t polydata<point_t>::insert(std::shared_ptr<geometry::geometric_object<point_t>> object)
        {
#ifdef __tpf_sanity_checks
            if (!(this->point_data.empty() && this->cell_data.empty() && this->object_data.empty()))
            {
                throw std::runtime_error(__tpf_error_message("Insertion of sole geometric objects not allowed if arrays are already attached."));
            }
#endif

            this->geometry.push_back(object);

            return this->geometry.size() - 1;
        }

        template <typename point_t>
        inline std::size_t polydata<point_t>::insert(std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& objects)
        {
#ifdef __tpf_sanity_checks
            if (!(this->point_data.empty() && this->cell_data.empty() && this->object_data.empty()))
            {
                throw std::runtime_error(__tpf_error_message("Insertion of sole geometric objects not allowed if arrays are already attached."));
            }
#endif

            this->geometry.insert(this->geometry.end(), objects.begin(), objects.end());

            return this->geometry.size() - 1;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
        inline std::size_t polydata<point_t>::insert(std::shared_ptr<geometry::geometric_object<point_t>> object,
            const polydata_element<value_t, rows, columns>& first_element, const data_elements_t&... data_elements)
        {
            const std::size_t num_elements = insert(object->get_num_points(), object->get_num_cells(), first_element, data_elements...);

#ifdef __tpf_sanity_checks
            if (num_elements != this->point_data.size() + this->cell_data.size() + this->object_data.size())
            {
                throw std::runtime_error(__tpf_error_message("Insertion of geometric failed. The number of data elements does not match the number of attached arrays."));
            }
#endif

            this->geometry.push_back(object);

            return this->geometry.size() - 1;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
        inline std::size_t polydata<point_t>::insert(const std::size_t num_points, const std::size_t num_cells,
            const polydata_element<value_t, rows, columns>& first_element, const data_elements_t&... data_elements)
        {
            std::size_t num_elements = 0;

            num_elements += insert(num_points, num_cells, first_element);
            num_elements += insert(num_points, num_cells, data_elements...);

            return num_elements;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline std::size_t polydata<point_t>::insert(const std::size_t num_points, const std::size_t num_cells, const polydata_element<value_t, rows, columns>& element)
        {
#ifdef __tpf_sanity_checks
            if ((element.topology == data::topology_t::POINT_DATA || element.topology == data::topology_t::TEXTURE_COORDINATES) && element.values.size() != num_points)
            {
                throw std::runtime_error(__tpf_error_message("Number of entries of the element doesn't match the number of points of the object."));
            }
            if (element.topology == data::topology_t::CELL_DATA && element.values.size() != num_cells)
            {
                throw std::runtime_error(__tpf_error_message("Number of entries of the element doesn't match the number of cells of the object."));
            }
            if (element.topology == data::topology_t::OBJECT_DATA && element.values.size() != 1)
            {
                throw std::runtime_error(__tpf_error_message("Number of entries of the element is not 1."));
            }
#endif

            if (element.topology == data::topology_t::POINT_DATA || element.topology == data::topology_t::TEXTURE_COORDINATES)
            {
#ifdef __tpf_sanity_checks
                if (!has_point_data(element.name))
                {
                    throw std::runtime_error(__tpf_error_message("Name not found in attached point data arrays."));
                }
#endif

                get_point_data_as<value_t, rows, columns>(element.name)->insert_back(element.values.begin(), element.values.end());
            }
            else if (element.topology == data::topology_t::CELL_DATA)
            {
#ifdef __tpf_sanity_checks
                if (!has_cell_data(element.name))
                {
                    throw std::runtime_error(__tpf_error_message("Name not found in attached cell data arrays."));
                }
#endif

                get_cell_data_as<value_t, rows, columns>(element.name)->insert_back(element.values.begin(), element.values.end());
            }
            else if (element.topology == data::topology_t::OBJECT_DATA)
            {
#ifdef __tpf_sanity_checks
                if (!has_object_data(element.name))
                {
                    throw std::runtime_error(__tpf_error_message("Name not found in attached object data arrays."));
                }
#endif
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }

            return 1;
        }

        template <typename point_t>
        inline void polydata<point_t>::add(std::shared_ptr<array_base> array, const topology_t type)
        {
            if (type == topology_t::POINT_DATA || type == topology_t::TEXTURE_COORDINATES)
            {
                if (array->get_num_elements() != get_num_points())
                {
                    throw std::runtime_error(__tpf_error_message("Number of entries in the array doesn't match the number of points."));
                }

                this->point_data.push_back(array);
            }
            else if (type == topology_t::CELL_DATA)
            {
                if (array->get_num_elements() != get_num_cells())
                {
                    throw std::runtime_error(__tpf_error_message("Number of entries in the array doesn't match the number of cells."));
                }

                this->cell_data.push_back(array);
            }
            else if (type == topology_t::OBJECT_DATA)
            {
                if (array->get_num_elements() != get_num_objects())
                {
                    throw std::runtime_error(__tpf_error_message("Number of entries in the array doesn't match the number of objects."));
                }

                this->object_data.push_back(array);
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }
        }

        template <typename point_t>
        void polydata<point_t>::merge(const polydata& other)
        {
            // Check that both poly data objects have the same arrays
            if (this->geometry.size() != 0)
            {
                if (this->point_data.size() != other.get_point_data().size())
                {
                    throw std::runtime_error(__tpf_error_message("Number of point data arrays does not match."));
                }
                if (this->cell_data.size() != other.get_cell_data().size())
                {
                    throw std::runtime_error(__tpf_error_message("Number of cell data arrays does not match."));
                }
                if (this->object_data.size() != other.get_object_data().size())
                {
                    throw std::runtime_error(__tpf_error_message("Number of object data arrays does not match."));
                }

                if (this->point_data.size() > 0)
                {
                    std::map<std::string, std::size_t> occurance;

                    for (const auto arr : this->point_data)
                    {
                        occurance[arr->get_name()] = 1;
                    }

                    for (const auto arr : other.get_point_data())
                    {
                        if (occurance.find(arr->get_name()) == occurance.end())
                        {
                            occurance[arr->get_name()] = 1;
                        }
                        else
                        {
                            ++occurance[arr->get_name()];
                        }
                    }

                    for (const auto& kv : occurance)
                    {
                        if (kv.second != 2)
                        {
                            throw std::runtime_error(__tpf_error_message("Point data arrays do not match."));
                        }
                    }
                }
                if (this->cell_data.size() > 0)
                {
                    std::map<std::string, std::size_t> occurance;

                    for (const auto arr : this->cell_data)
                    {
                        occurance[arr->get_name()] = 1;
                    }

                    for (const auto arr : other.get_cell_data())
                    {
                        if (occurance.find(arr->get_name()) == occurance.end())
                        {
                            occurance[arr->get_name()] = 1;
                        }
                        else
                        {
                            ++occurance[arr->get_name()];
                        }
                    }

                    for (const auto& kv : occurance)
                    {
                        if (kv.second != 2)
                        {
                            throw std::runtime_error(__tpf_error_message("Cell data arrays do not match."));
                        }
                    }
                }
                if (this->object_data.size() > 0)
                {
                    std::map<std::string, std::size_t> occurance;

                    for (const auto arr : this->object_data)
                    {
                        occurance[arr->get_name()] = 1;
                    }

                    for (const auto arr : other.get_object_data())
                    {
                        if (occurance.find(arr->get_name()) == occurance.end())
                        {
                            occurance[arr->get_name()] = 1;
                        }
                        else
                        {
                            ++occurance[arr->get_name()];
                        }
                    }

                    for (const auto& kv : occurance)
                    {
                        if (kv.second != 2)
                        {
                            throw std::runtime_error(__tpf_error_message("Object data arrays do not match."));
                        }
                    }
                }
            }

            // Copy geometry
            for (const auto object : other.get_geometry())
            {
                this->geometry.push_back(object->clone());
            }

            // Copy arrays
            for (auto arr : this->point_data)
            {
                arr->merge(other.get_point_data(arr->get_name()));
            }
            for (auto arr : this->cell_data)
            {
                arr->merge(other.get_cell_data(arr->get_name()));
            }
            for (auto arr : this->object_data)
            {
                arr->merge(other.get_object_data(arr->get_name()));
            }
        }

        template <typename point_t>
        inline std::shared_ptr<polydata<point_t>> polydata<point_t>::clone(const math::transformer<point_t, 3>& trafo) const
        {
            auto copy = std::make_shared<polydata<point_t>>(*this);

            for (auto& object : copy->get_geometry())
            {
                object->transform(trafo);
            }

            return copy;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline std::shared_ptr<array<value_t, rows, columns>> polydata<point_t>::create(const std::string& name, const topology_t type)
        {
            std::shared_ptr<array<value_t, rows, columns>> new_array = nullptr;
            if (type == topology_t::POINT_DATA || type == topology_t::TEXTURE_COORDINATES)
            {
                new_array = std::make_shared<array<value_t, rows, columns>>(name, get_num_points());
            }
            else if (type == topology_t::CELL_DATA)
            {
                new_array = std::make_shared<array<value_t, rows, columns>>(name, get_num_cells());
            }
            else if (type == topology_t::OBJECT_DATA)
            {
                new_array = std::make_shared<array<value_t, rows, columns>>(name, get_num_objects());
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }

            add(new_array, type);

            return new_array;
        }

        template <typename point_t>
        inline polydata<point_t> polydata<point_t>::extract_area(const area_t<point_t>& subarea, const bool include_arrays) const
        {
            const auto indices = extract_indices(subarea);

            auto geometry = extract_geometry(indices);

            polydata<point_t> extraction(geometry);

            if (include_arrays)
            {
                for (const auto& data_array : this->point_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::POINT_DATA);
                }
                for (const auto& data_array : this->cell_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::CELL_DATA);
                }
                for (const auto& data_array : this->object_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::OBJECT_DATA);
                }
            }

            return extraction;
        }

        template <typename point_t>
        inline polydata<point_t> polydata<point_t>::extract_area(const area_t<point_t>& subarea, const bool include_arrays)
        {
            const auto indices = extract_indices(subarea);

            auto geometry = extract_geometry(indices);

            polydata<point_t> extraction(geometry);

            if (include_arrays)
            {
                for (const auto& data_array : this->point_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::POINT_DATA);
                }
                for (const auto& data_array : this->cell_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::CELL_DATA);
                }
                for (const auto& data_array : this->object_data)
                {
                    extraction.add(data_array->extract(indices), topology_t::OBJECT_DATA);
                }
            }

            return extraction;
        }

        template <typename point_t>
        template <typename... data_elements_t>
        inline polydata<point_t> polydata<point_t>::extract_area(const area_t<point_t>& subarea, const data_elements_t&... data_elements) const
        {
            const auto indices = extract_indices(subarea);

            auto geometry = extract_geometry(indices);
            auto arrays = extract_array(data_elements...);

            polydata<point_t> extraction(geometry);

            for (auto data_array : arrays)
            {
                extraction.add(data_array.first, data_array.second);
            }

            return extraction;
        }

        template <typename point_t>
        template <typename... data_elements_t>
        inline polydata<point_t> polydata<point_t>::extract_area(const area_t<point_t>& subarea, const data_elements_t&... data_elements)
        {
            const auto indices = extract_indices(subarea);

            auto geometry = extract_geometry(indices);
            auto arrays = extract_array(data_elements...);

            polydata<point_t> extraction(geometry);

            for (auto data_array : arrays)
            {
                extraction.add(data_array.first, data_array.second);
            }

            return extraction;
        }

        template <typename point_t>
        inline std::vector<std::size_t> polydata<point_t>::extract_indices(const area_t<point_t>& subarea) const
        {
            std::vector<std::size_t> indices;
            indices.reserve(this->geometry.size());

            for (std::size_t i = 0; i < this->geometry.size(); ++i)
            {
                const auto points = this->geometry[i]->get_points();
                bool inside = false;

                for (const auto& point : points)
                {
                    inside |=
                        subarea[0].first <= point[0] && point[0] <= subarea[0].second &&
                        subarea[1].first <= point[1] && point[1] <= subarea[1].second &&
                        subarea[2].first <= point[2] && point[2] <= subarea[2].second;
                }

                if (inside)
                {
                    indices.push_back(i);
                }
            }

            return indices;
        }

        template <typename point_t>
        inline std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> polydata<point_t>::extract_geometry(const std::vector<std::size_t>& indices) const
        {
            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> extraction(indices.size());

            for (std::size_t i = 0; i < indices.size(); ++i)
            {
                extraction[i] = this->geometry[indices[i]]->clone();
            }

            return extraction;
        }

        template <typename point_t>
        inline std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> polydata<point_t>::extract_geometry(const std::vector<std::size_t>& indices)
        {
            std::vector<std::shared_ptr<geometry::geometric_object<point_t>>> extraction(indices.size());

            for (std::size_t i = 0; i < indices.size(); ++i)
            {
                extraction[i] = this->geometry[indices[i]];
            }

            return extraction;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns, typename... data_elements_t>
        inline std::vector<std::pair<std::shared_ptr<data::array_base>, data::topology_t>> polydata<point_t>::extract_array(
            const std::vector<std::size_t>& indices, const data_information<value_t, rows, columns>& first_element, const data_elements_t&... data_elements) const
        {
            auto first_array = extract_array(first_element);
            auto further_arrays = extract_array(data_elements...);

            further_arrays.insert(further_arrays.end(), first_array.begin(), first_array.end());

            return further_arrays;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline std::vector<std::pair<std::shared_ptr<data::array_base>, data::topology_t>> polydata<point_t>::extract_array(
            const std::vector<std::size_t>& indices, const data_information<value_t, rows, columns>& element) const
        {
            std::vector<std::pair<std::shared_ptr<data::array_base>, data::topology_t>> extraction;

            if (element.topology == topology_t::POINT_DATA || element.topology == topology_t::TEXTURE_COORDINATES)
            {
                auto extracted_array = get_point_data(element.name)->extract(indices);
                extraction.push_back(std::make_pair(extracted_array, topology_t::POINT_DATA));
            }
            else if (element.topology == topology_t::CELL_DATA)
            {
                auto extracted_array = get_cell_data(element.name)->extract(indices);
                extraction.push_back(std::make_pair(extracted_array, topology_t::CELL_DATA));
            }
            else if (element.topology == topology_t::OBJECT_DATA)
            {
                auto extracted_array = get_object_data(element.name)->extract(indices);
                extraction.push_back(std::make_pair(extracted_array, topology_t::OBJECT_DATA));
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }

            return extraction;
        }

        template <typename point_t>
        inline const std::shared_ptr<geometry::geometric_object<point_t>> polydata<point_t>::get_object(const std::size_t id) const
        {
#ifdef __tpf_range_checks
            if (id >= this->geometry.size())
            {
                throw std::runtime_error(__tpf_error_message("Geometric object with given ID could not be found."));
            }
#endif

            return this->geometry[id];
        }

        template <typename point_t>
        inline const std::vector<std::shared_ptr<geometry::geometric_object<point_t>>>& polydata<point_t>::get_geometry() const
        {
            return this->geometry;
        }

        template <typename point_t>
        bool polydata<point_t>::has_point_data(const std::string& name) const
        {
            for (auto entry : this->point_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    return true;
                }
            }

            return false;
        }

        template <typename point_t>
        inline const std::vector<std::shared_ptr<array_base>>& polydata<point_t>::get_point_data() const
        {
            return this->point_data;
        }

        template <typename point_t>
        inline const std::shared_ptr<array_base> polydata<point_t>::get_point_data(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->point_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find point data '", name, "'."));
            }

            return ret_value;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline const std::shared_ptr<array<value_t, rows, columns>> polydata<point_t>::get_point_data_as(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->point_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find point data '", name, "'."));
            }

            return std::dynamic_pointer_cast<array<value_t, rows, columns>, array_base>(ret_value);
        }

        template <typename point_t>
        bool polydata<point_t>::has_cell_data(const std::string& name) const
        {
            for (auto entry : this->cell_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    return true;
                }
            }

            return false;
        }

        template <typename point_t>
        inline const std::vector<std::shared_ptr<array_base>>& polydata<point_t>::get_cell_data() const
        {
            return this->cell_data;
        }

        template <typename point_t>
        inline const std::shared_ptr<array_base> polydata<point_t>::get_cell_data(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->cell_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find cell data '", name, "'."));
            }

            return ret_value;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline const std::shared_ptr<array<value_t, rows, columns>> polydata<point_t>::get_cell_data_as(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->cell_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find cell data '", name, "'."));
            }

            return std::dynamic_pointer_cast<array<value_t, rows, columns>, array_base>(ret_value);
        }

        template <typename point_t>
        bool polydata<point_t>::has_object_data(const std::string& name) const
        {
            for (auto entry : this->object_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    return true;
                }
            }

            return false;
        }

        template <typename point_t>
        inline const std::vector<std::shared_ptr<array_base>>& polydata<point_t>::get_object_data() const
        {
            return this->object_data;
        }

        template <typename point_t>
        inline const std::shared_ptr<array_base> polydata<point_t>::get_object_data(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->object_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find object data '", name, "'."));
            }

            return ret_value;
        }

        template <typename point_t>
        template <typename value_t, std::size_t rows, std::size_t columns>
        inline const std::shared_ptr<array<value_t, rows, columns>> polydata<point_t>::get_object_data_as(const std::string& name) const
        {
            std::shared_ptr<array_base> ret_value = nullptr;

            for (auto entry : this->object_data)
            {
                if (entry->get_name().compare(name) == 0)
                {
                    ret_value = entry;
                }
            }

            if (ret_value == nullptr)
            {
                throw std::runtime_error(__tpf_error_message("Could not find object data '", name, "'."));
            }

            return std::dynamic_pointer_cast<array<value_t, rows, columns>, array_base>(ret_value);
        }

        template <typename point_t>
        inline std::size_t polydata<point_t>::get_num_points() const
        {
            std::size_t num_points = 0;

            for (const auto obj : this->geometry)
            {
                num_points += obj->get_num_points();
            }

            return num_points;
        }

        template <typename point_t>
        inline std::size_t polydata<point_t>::get_num_cells() const
        {
            std::size_t num_cells = 0;

            for (const auto obj : this->geometry)
            {
                num_cells += obj->get_num_cells();
            }

            return num_cells;
        }

        template <typename point_t>
        inline std::size_t polydata<point_t>::get_num_objects() const
        {
            return this->geometry.size();
        }

        template <typename point_t>
        inline std::vector<std::size_t> polydata<point_t>::get_num_object_cells() const
        {
            std::vector<std::size_t> num_cells(this->geometry.size());
            for (std::size_t i = 0; i < this->geometry.size(); ++i)
            {
                num_cells[i] = this->geometry[i]->get_num_cells();
            }
            return num_cells;
        }
    }
}
