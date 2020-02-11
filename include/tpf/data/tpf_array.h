#pragma once

#include "tpf_array_base.h"

#include "Eigen/Dense"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Array storing arbitrary types or mathematical vectors and matrices
        /// </summary>
        /// <template name="value_t">Value type of the data or of the vector/matrix</template>
        /// <template name="rows">Number of rows of the vector/matrix</template>
        /// <template name="columns">Number of columns of the vector/matrix</template>
        template <typename value_t, int rows = 1, int columns = 1>
        class array : public array_base
        {
        public:
            static_assert(rows > 0, "Number of rows must be larger than zero");
            static_assert(columns > 0, "Number of columns must be larger than zero");

            /// Used value types
            using value_type = value_t;
            using element_type = typename std::conditional<(rows > 1 || columns > 1), Eigen::Matrix<value_t, rows, columns>, value_t>::type;
            using return_type = typename std::conditional<(rows > 1 || columns > 1), Eigen::Map<element_type>, element_type&>::type;

            /// <summary>
            /// Create array with given number of elements
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="num_elements">Number of elements</param>
            array(std::string name, std::size_t num_elements = 0);

            /// <summary>
            /// Create array with given data
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="data">Data</param>
            array(std::string name, std::vector<value_t> data);

            /// <summary>
            /// Default copy/move constructor/assignment
            /// </summary>
            array(const array&) = default;
            array(array&&) noexcept = default;
            array& operator=(const array&) = default;
            array& operator=(array&&) noexcept = default;

            /// <summary>
            /// Clone object (deep copy)
            /// </summary>
            /// <returns>Cloned object</returns>
            virtual std::shared_ptr<array_base> clone() const override;

            /// <summary>
            /// Merge other array into this one, by appending its elements
            /// </summary>
            /// <param name="other">Other array to copy from</param>
            virtual void merge(const std::shared_ptr<array_base> other) override;

            /// <summary>
            /// Extract specific elements from this array
            /// </summary>
            /// <param name="indices">Indices of elements to extract</param>
            /// <returns>Extracted array with specified elements</returns>
            virtual std::shared_ptr<array_base> extract(const std::vector<std::size_t>& indices) const override;

            /// <summary>
            /// Initialize array with given value
            /// </summary>
            /// <param name="value">Value for every entry</param>
            void initialize(const value_t& value) noexcept;

            /// <summary>
            /// Initialize array with given vector
            /// </summary>
            /// <param name="value">Vector or matrix for every entry</param>
            template <int _enable_rows = rows, int _enable_columns = columns>
            void initialize(const element_type& value, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type* = nullptr) noexcept;

            /// <summary>
            /// Get number of components
            /// </summary>
            /// <return>Number of components</return>
            constexpr std::size_t get_num_components() const noexcept;

            /// <summary>
            /// Get number of rows
            /// </summary>
            /// <return>Number of rows</return>
            constexpr std::size_t get_num_rows() const noexcept;

            /// <summary>
            /// Get number of columns
            /// </summary>
            /// <return>Number of columns</return>
            constexpr std::size_t get_num_columns() const noexcept;

            /// <summary>
            /// Get array size
            /// </summary>
            /// <return>Array size</return>
            std::size_t get_size() const noexcept;

            /// <summary>
            /// Access data
            /// </summary>
            /// <return>Data</return>
            const std::vector<value_t>& get_data() const noexcept;

            /// <summary>
            /// Access data
            /// </summary>
            /// <return>Data</return>
            std::vector<value_t>& get_data() noexcept;

            /// <summary>
            /// Set data
            /// </summary>
            /// <param name="data">Data</param>
            void set_data(std::vector<value_t> data);

            /// <summary>
            /// Access data component at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <param name="component">Component index</param>
            /// <return>Component of an entry</return>
            const value_t& operator()(std::size_t element, std::size_t component) const;

            /// <summary>
            /// Access data component at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <param name="component">Component index</param>
            /// <return>Component of an entry</return>
            value_t& operator()(std::size_t element, std::size_t component);

            /// <summary>
            /// Access data component at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <param name="component">Component index</param>
            /// <return>Component of an entry</return>
            const value_t& at(std::size_t element, std::size_t component) const;

            /// <summary>
            /// Access data component at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <param name="component">Component index</param>
            /// <return>Component of an entry</return>
            value_t& at(std::size_t element, std::size_t component);

            /// <summary>
            /// Access data at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <return>Vector or scalar value</return>
            const return_type operator()(std::size_t element) const;

            /// <summary>
            /// Access data at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <return>Vector or scalar value</return>
            return_type operator()(std::size_t element);

            /// <summary>
            /// Access data at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <return>Vector or scalar value</return>
            const return_type at(std::size_t element) const;

            /// <summary>
            /// Access data at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <return>Vector or scalar value</return>
            return_type at(std::size_t element);

            /// <summary>
            /// Extend array by adding one element at the end
            /// </summary>
            /// <param name="new_element">New element</param>
            template <int _enable_rows = rows, int _enable_columns = columns>
            void push_back(element_type new_element, typename std::enable_if<(_enable_rows == 1 && _enable_columns == 1)>::type* = nullptr);

            template <int _enable_rows = rows, int _enable_columns = columns>
            void push_back(const element_type& new_element, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type* = nullptr);

            /// <summary>
            /// Extend array by adding elements at the end
            /// </summary>
            /// <param name="new_entries">New entries</param>
            template <typename iterator>
            void insert_back(iterator first, iterator end);

            /// <summary>
            /// Reserve space for given number of elements (does not shrink the array)
            /// </summary>
            /// <param name="num_elements">Number of elements reserved</param>
            void reserve(std::size_t num_elements);

            /// <summary>
            /// Resize data for given number of elements (may shrink the array)
            /// </summary>
            /// <param name="num_elements">New number of elements</param>
            void resize(std::size_t num_elements);

        private:
            /// <summary>
            /// Access data at the given element index
            /// </summary>
            /// <param name="element">Element index</param>
            /// <return>Vector or scalar value</return>
            template <int _enable_rows = rows, int _enable_columns = columns>
            return_type access(std::size_t element, typename std::enable_if<(_enable_rows == 1 && _enable_columns == 1)>::type* = nullptr);

            template <int _enable_rows = rows, int _enable_columns = columns>
            return_type access(std::size_t element, typename std::enable_if<(_enable_rows > 1 || _enable_columns > 1)>::type* = nullptr);

            /// Data
            std::vector<value_t> data;
        };
    }
}

#include "tpf_array.inl"