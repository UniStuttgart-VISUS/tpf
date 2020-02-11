#include "tpf_data.h"

#include "tpf_vtk_traits.h"

#include "../data/tpf_array.h"
#include "../data/tpf_data_information.h"

#include "../log/tpf_log.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace vtk
    {
        template <typename expected_t, typename found_t>
        inline void warn_mismatch(const std::string& array_name)
        {
#ifdef __tpf_debug
            if (sizeof(expected_t) < sizeof(found_t))
            {
                log::warning_message(__tpf_warning_message("Warning: Data type mismatch in array '", array_name, "'! ",
                    typeid(expected_t).name(), " expected, ", typeid(found_t).name(), " found.\n",
                    "Possible narrowing conversion!"));
            }
            else
            {
                log::warning_message(__tpf_warning_message("Warning: Data type mismatch in array '", array_name, "'! ",
                    typeid(expected_t).name(), " expected, ", typeid(found_t).name(), " found."));
            }
#endif
        }

        inline bool has_data(vtkDataSet* data, const data::topology_t data_type)
        {
            if (data_type == data::topology_t::CELL_DATA)
            {
                return data->GetCellData()->GetNumberOfArrays() > 0;
            }
            else if (data_type == data::topology_t::POINT_DATA)
            {
                return data->GetPointData()->GetNumberOfArrays() > 0;
            }
            else if (data_type == data::topology_t::OBJECT_DATA)
            {
                // Object data topology is not supported on vtk data.
                return false;
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }
        }

        inline bool has_array(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name)
        {
            if (data_type == data::topology_t::CELL_DATA)
            {
                return data->GetCellData()->GetArray(array_name.c_str()) != nullptr;
            }
            else if (data_type == data::topology_t::POINT_DATA)
            {
                return data->GetPointData()->GetArray(array_name.c_str()) != nullptr;
            }
            else if (data_type == data::topology_t::OBJECT_DATA)
            {
                // Object data topology is not supported on vtk data.
                return false;
            }
            else
            {
                throw std::runtime_error(__tpf_error_message("Invalid topology type."));
            }
        }

        template <typename value_t, class array_t, bool dynamic>
        inline std::vector<value_t> get_data(vtkAbstractArray* data, const std::string& array_name)
        {
            try
            {
                static_assert(std::is_base_of<vtkDataArray, array_t>::value, "Class must be a derived class from vtkDataArray!");

                // Try downcasting to the given data type
                auto temp = array_t::SafeDownCast(data);

                if (temp != nullptr)
                {
                    std::vector<value_t> data_array;
                    data_array.reserve(static_cast<std::size_t>(data->GetNumberOfValues()));

                    for (int i = 0; i < data->GetNumberOfValues(); ++i)
                    {
                        data_array.push_back(static_cast<value_t>(temp->GetValue(i)));
                    }

                    return data_array;
                }
                else if (dynamic)
                {
                    // Downcast to different data types but issue a warning
                    if (std::is_floating_point<value_t>::value)
                    {
                        // Try double and float
                        switch (data->GetDataType())
                        {
                        case VTK_FLOAT:
                            warn_mismatch<value_t, float>(array_name);

                            return get_data<value_t, typename vtk_array<float>::type, false>(data, array_name);
                        case VTK_DOUBLE:
                            warn_mismatch<value_t, double>(array_name);

                            return get_data<value_t, typename vtk_array<double>::type, false>(data, array_name);
                        }
                    }
                    else if (std::is_integral<value_t>::value && std::is_signed<value_t>::value)
                    {
                        // Try signed integer types
                        switch (data->GetDataType())
                        {
                        case VTK_CHAR:
                            warn_mismatch<value_t, char>(array_name);

                            return get_data<value_t, typename vtk_array<char>::type, false>(data, array_name);
                        case VTK_INT:
                            warn_mismatch<value_t, int>(array_name);

                            return get_data<value_t, typename vtk_array<int>::type, false>(data, array_name);
                        case VTK_LONG:
                            warn_mismatch<value_t, long>(array_name);

                            return get_data<value_t, typename vtk_array<long>::type, false>(data, array_name);
                        case VTK_LONG_LONG:
                            warn_mismatch<value_t, long long>(array_name);

                            return get_data<value_t, typename vtk_array<long long>::type, false>(data, array_name);
                        case VTK_SHORT:
                            warn_mismatch<value_t, short>(array_name);

                            return get_data<value_t, typename vtk_array<short>::type, false>(data, array_name);
                        case VTK_SIGNED_CHAR:
                            warn_mismatch<value_t, signed char>(array_name);

                            return get_data<value_t, typename vtk_array<signed char>::type, false>(data, array_name);
                        }
                    }
                    else if (std::is_integral<value_t>::value && std::is_unsigned<value_t>::value)
                    {
                        // Try unsigned integer types
                        switch (data->GetDataType())
                        {
                        case VTK_UNSIGNED_CHAR:
                            warn_mismatch<value_t, unsigned char>(array_name);

                            return get_data<value_t, typename vtk_array<unsigned char>::type, false>(data, array_name);
                        case VTK_UNSIGNED_INT:
                            warn_mismatch<value_t, unsigned int>(array_name);

                            return get_data<value_t, typename vtk_array<unsigned int>::type, false>(data, array_name);
                        case VTK_UNSIGNED_LONG:
                            warn_mismatch<value_t, unsigned long>(array_name);

                            return get_data<value_t, typename vtk_array<unsigned long>::type, false>(data, array_name);
                        case VTK_UNSIGNED_LONG_LONG:
                            warn_mismatch<value_t, unsigned long long>(array_name);

                            return get_data<value_t, typename vtk_array<unsigned long long>::type, false>(data, array_name);
                        case VTK_UNSIGNED_SHORT:
                            warn_mismatch<value_t, unsigned short>(array_name);

                            return get_data<value_t, typename vtk_array<unsigned short>::type, false>(data, array_name);
                        }
                    }

                    // Throw exception for incompatible types
                    throw std::runtime_error(__tpf_error_message("Given type and found VTK array type are incompatible."));
                }
                else
                {
                    // Throw exception for incompatible types
                    throw std::runtime_error(__tpf_error_message("Given type and found VTK array type are incompatible."));
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Failed getting data."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Failed getting data."));
            }
        }

        template <typename value_t, class array_t, bool dynamic>
        inline std::vector<value_t> get_data(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name)
        {
            try
            {
                static_assert(std::is_base_of<vtkDataArray, array_t>::value, "Class must be a derived class from vtkDataArray!");

                if (data_type == data::topology_t::CELL_DATA)
                {
                    return get_data<value_t, array_t, dynamic>(data->GetCellData()->GetAbstractArray(array_name.c_str()), array_name.c_str());
                }
                else if (data_type == data::topology_t::POINT_DATA)
                {
                    return get_data<value_t, array_t, dynamic>(data->GetPointData()->GetAbstractArray(array_name.c_str()), array_name.c_str());
                }
                else if (data_type == data::topology_t::OBJECT_DATA)
                {
                    throw std::runtime_error(__tpf_error_message("Object data topology is not supported on vtk data."));
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Invalid topology type."));
                }
            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Failed getting data."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Failed getting data."));
            }
        }

        template <typename value_t, class array_t>
        inline void set_data(vtkAbstractArray* data, const std::vector<value_t>& data_array, const std::size_t num_components)
        {
            auto temp_array = array_t::SafeDownCast(data);
            temp_array->SetNumberOfComponents(static_cast<int>(num_components));
            temp_array->SetNumberOfTuples(static_cast<int>(data_array.size() / num_components));

            std::transform(data_array.begin(), data_array.end(), CHECKED_ITERATOR(temp_array, array_t), [](const value_t& item) { return item; });
        }

        template <typename value_t, class array_t>
        inline void set_data(vtkDataSet* data, const data::topology_t data_type, const std::string& array_name,
            const std::vector<value_t>& data_array, const std::size_t num_components)
        {
            try
            {
                static_assert(std::is_base_of<vtkDataArray, array_t>::value, "Class must be a derived class from vtkDataArray!");

                auto temp_array = vtkSmartPointer<array_t>::New();
                temp_array->SetName(array_name.c_str());

                set_data<value_t, array_t>(temp_array, data_array, num_components);

                if (data_type == data::topology_t::CELL_DATA)
                {
                    data->GetCellData()->AddArray(temp_array);
                }
                else if (data_type == data::topology_t::POINT_DATA)
                {
                    data->GetPointData()->AddArray(temp_array);
                }
                else if (data_type == data::topology_t::OBJECT_DATA)
                {
                    throw std::runtime_error(__tpf_error_message("Object data topology is not supported on vtk data."));
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Invalid topology type."));
                }

            }
            catch (const std::runtime_error& ex)
            {
                throw std::runtime_error(__tpf_nested_error_message(ex.what(), "Failed setting data."));
            }
            catch (...)
            {
                throw std::runtime_error(__tpf_error_message("Failed setting data."));
            }
        }

        template <typename value_t, class array_t, int rows, int columns>
        inline void set_data(vtkDataSet* data, const data::topology_t data_type, const data::array<value_t, rows, columns>& array)
        {
            set_data(data, data_type, array.get_name(), array.get_data(), array.get_num_components());
        }
    }
}
