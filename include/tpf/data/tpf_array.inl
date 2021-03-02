#include "tpf_array.h"

#include "tpf_array_base.h"

#include "../log/tpf_log.h"

#include "Eigen/Dense"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace tpf
{
    namespace data
    {
        template <typename value_t, int rows, int columns>
        inline array<value_t, rows, columns>::array(std::string name, const std::size_t num_elements) : array_base(std::move(name), num_elements)
        {
            this->data.resize(num_elements * rows * columns);
        }

        template <typename value_t, int rows, int columns>
        inline array<value_t, rows, columns>::array(std::string name, std::vector<value_t> data) : array_base(std::move(name), data.size() / (rows * columns))
        {
            if (data.size() % (rows * columns) != static_cast<std::size_t>(0))
            {
                throw std::runtime_error(__tpf_error_message("Data size doesn't match the number of components."));
            }

            std::swap(this->data, data);
        }

        template <typename value_t, int rows, int columns>
        inline std::shared_ptr<array_base> array<value_t, rows, columns>::clone() const
        {
            return std::make_shared<array<value_t, rows, columns>>(*this);
        }

        template <typename value_t, int rows, int columns>
        inline void array<value_t, rows, columns>::merge(const std::shared_ptr<array_base> other)
        {
            const auto other_array = std::dynamic_pointer_cast<array<value_t, rows, columns>>(other);

            this->data.insert(this->data.end(), other_array->get_data().begin(), other_array->get_data().end());
        }

        template <typename value_t, int rows, int columns>
        inline std::shared_ptr<array_base> array<value_t, rows, columns>::extract(const std::vector<std::size_t>& indices) const
        {
            auto extraction = std::make_shared<array<value_t, rows, columns>>(array_base::name, indices.size());

            for (std::size_t i = 0; i < indices.size(); ++i)
            {
                (*extraction)(i) = operator()(indices[i]);
            }

            return std::dynamic_pointer_cast<array_base>(extraction);
        }

        template <typename value_t, int rows, int columns>
        inline void array<value_t, rows, columns>::initialize(const value_t& value) noexcept
        {
            std::fill(this->data.begin(), this->data.end(), value);
        }

        template <typename value_t, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline void array<value_t, rows, columns>::initialize(const element_type& value, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type*) noexcept
        {
            for (std::size_t index = 0; index < this->data.size(); index += rows * columns)
            {
                Eigen::Map<element_type> mapped(&this->data[index]);
                mapped = value;
            }
        }

        template <typename value_t, int rows, int columns>
        inline std::size_t array<value_t, rows, columns>::get_num_components_dynamic() const noexcept
        {
            return rows * columns;
        }

        template <typename value_t, int rows, int columns>
        inline constexpr std::size_t array<value_t, rows, columns>::get_num_components() const noexcept
        {
            return rows * columns;
        }

        template <typename value_t, int rows, int columns>
        inline constexpr std::size_t array<value_t, rows, columns>::get_num_rows() const noexcept
        {
            return rows;
        }

        template <typename value_t, int rows, int columns>
        inline constexpr std::size_t array<value_t, rows, columns>::get_num_columns() const noexcept
        {
            return columns;
        }

        template <typename value_t, int rows, int columns>
        inline std::size_t array<value_t, rows, columns>::get_size() const noexcept
        {
            return this->data.size();
        }

        template <typename value_t, int rows, int columns>
        inline const std::vector<double>& array<value_t, rows, columns>::get_data_dynamic() const
        {
            return get_data_dynamic_impl();
        }

        template <typename value_t, int rows, int columns>
        template <typename local_value_t>
        const std::vector<double>& array<value_t, rows, columns>::get_data_dynamic_impl(
            typename std::enable_if<std::is_same<double, local_value_t>::value>::type*) const
        {
            return this->data;
        }

        template <typename value_t, int rows, int columns>
        template <typename local_value_t>
        const std::vector<double>& array<value_t, rows, columns>::get_data_dynamic_impl(
            typename std::enable_if<!std::is_same<double, local_value_t>::value>::type*) const
        {
            throw std::runtime_error(__tpf_error_message("Dynamic data access only available for double arrays."));
        }

        template <typename value_t, int rows, int columns>
        inline const std::vector<value_t>& array<value_t, rows, columns>::get_data() const noexcept
        {
            return this->data;
        }

        template <typename value_t, int rows, int columns>
        inline std::vector<value_t>& array<value_t, rows, columns>::get_data() noexcept
        {
            return this->data;
        }

        template <typename value_t, int rows, int columns>
        inline void array<value_t, rows, columns>::set_data(std::vector<value_t> data)
        {
            if (data.size() != this->data.size())
            {
                throw std::runtime_error(__tpf_error_message("New data size must match the local data size."));
            }

            std::swap(this->data, data);
        }

        template <typename value_t, int rows, int columns>
        inline const value_t& array<value_t, rows, columns>::operator()(const std::size_t element, const std::size_t component) const
        {
#ifdef __tpf_range_checks
            if (element >= array_base::num_elements)
            {
                throw std::runtime_error(__tpf_error_message("Element index must be smaller than the number of elements available."));
            }

            if (component >= rows * columns)
            {
                throw std::runtime_error(__tpf_error_message("Component index must be smaller than the number of components available."));
            }
#endif

            return this->data[element * rows * columns + component];
        }

        template <typename value_t, int rows, int columns>
        inline value_t& array<value_t, rows, columns>::operator()(const std::size_t element, const std::size_t component)
        {
#ifdef __tpf_range_checks
            if (element >= array_base::num_elements)
            {
                throw std::runtime_error(__tpf_error_message("Element index must be smaller than the number of elements available."));
            }

            if (component >= rows * columns)
            {
                throw std::runtime_error(__tpf_error_message("Component index must be smaller than the number of components available."));
            }
#endif

            return this->data[element * rows * columns + component];
        }

        template <typename value_t, int rows, int columns>
        inline const value_t& array<value_t, rows, columns>::at(const std::size_t element, const std::size_t component) const
        {
            return this->operator()(element, component);
        }

        template <typename value_t, int rows, int columns>
        inline value_t& array<value_t, rows, columns>::at(const std::size_t element, const std::size_t component)
        {
            return this->operator()(element, component);
        }

        template <typename value_t, int rows, int columns>
        inline const typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::operator()(const std::size_t element) const
        {
#ifdef __tpf_range_checks
            if (element >= array_base::num_elements)
            {
                throw std::runtime_error(__tpf_error_message("Element index must be smaller than the number of elements available."));
            }
#endif

            return const_cast<array<value_t, rows, columns>*>(this)->access(element);
        }

        template <typename value_t, int rows, int columns>
        inline typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::operator()(const std::size_t element)
        {
#ifdef __tpf_range_checks
            if (element >= array_base::num_elements)
            {
                throw std::runtime_error(__tpf_error_message("Element index must be smaller than the number of elements available."));
            }
#endif

            return access(element);
        }

        template <typename value_t, int rows, int columns>
        inline const typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::at(const std::size_t element) const
        {
            return this->operator()(element);
        }

        template <typename value_t, int rows, int columns>
        inline typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::at(const std::size_t element)
        {
            return this->operator()(element);
        }

        template <typename value_t, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::access(std::size_t element,
            typename std::enable_if<(_enable_rows == 1 && _enable_columns == 1)>::type*)
        {
            return this->data[element];
        }

        template <typename value_t, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline typename array<value_t, rows, columns>::return_type array<value_t, rows, columns>::access(std::size_t element,
            typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type*)
        {
            return return_type(&this->data[element * rows * columns]);
        }

        template <typename value_t, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline void array<value_t, rows, columns>::push_back(element_type new_element, typename std::enable_if<(_enable_rows == 1 && _enable_columns == 1)>::type*)
        {
            this->data.push_back(std::move(new_element));

            ++array_base::num_elements;
        }

        template <typename value_t, int rows, int columns>
        template <int _enable_rows, int _enable_columns>
        inline void array<value_t, rows, columns>::push_back(const element_type& new_element, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type*)
        {
            for (std::size_t i = 0; i < rows * columns; ++i)
            {
                this->data.push_back(new_element(i));
            }

            ++array_base::num_elements;
        }

        template <typename value_t, int rows, int columns>
        template <typename iterator>
        inline void array<value_t, rows, columns>::insert_back(iterator current, const iterator end)
        {
            for (; current != end; ++current)
            {
                push_back(*current);
            }
        }

        template <typename value_t, int rows, int columns>
        inline void array<value_t, rows, columns>::reserve(const std::size_t num_elements)
        {
            this->data.reserve(num_elements * rows * columns);
        }

        template <typename value_t, int rows, int columns>
        inline void array<value_t, rows, columns>::resize(const std::size_t num_elements)
        {
            array_base::num_elements = num_elements;

            this->data.resize(num_elements * rows * columns);
        }
    }
}