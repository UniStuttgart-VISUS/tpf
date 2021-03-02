#pragma once

#include <memory>
#include <string>
#include <vector>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Base class for array data structure
        /// </summary>
        class array_base
        {
        public:
            /// <summary>
            /// Create a named array with initial number of elements
            /// </summary>
            /// <param name="name">Array name</param>
            /// <param name="num_elements">Number of initial elements</param>
            array_base(std::string name, std::size_t num_elements) noexcept;

            /// <summary>
            /// Default copy/move constructor/assignment
            /// </summary>
            array_base(const array_base&) = default;
            array_base(array_base&&) noexcept = default;
            array_base& operator=(const array_base&) = default;
            array_base& operator=(array_base&&) noexcept = default;

            /// <summary>
            /// Clone object (deep copy)
            /// </summary>
            /// <returns>Cloned object</returns>
            virtual std::shared_ptr<array_base> clone() const = 0;

            /// <summary>
            /// Merge other array into this one, by appending its elements
            /// </summary>
            /// <param name="other">Array to copy from</param>
            virtual void merge(const std::shared_ptr<array_base> other) = 0;

            /// <summary>
            /// Extract specific elements from this array
            /// </summary>
            /// <param name="indices">Indices of elements to extract</param>
            /// <returns>Extracted array with specified elements</returns>
            virtual std::shared_ptr<array_base> extract(const std::vector<std::size_t> &indices) const = 0;

            /// <summary>
            /// Get the array name
            /// </summary>
            /// <return>Array name</return>
            const std::string& get_name() const noexcept;

            /// <summary>
            /// Set a new array name
            /// </summary>
            /// <param name="name">New array name</param>
            void set_name(std::string name) noexcept;

            /// <summary>
            /// Get number of stored elements
            /// </summary>
            /// <return>Number of stored elements</return>
            std::size_t get_num_elements() const noexcept;

            /// <summary>
            /// Get number of components
            /// </summary>
            /// <return>Number of stored elements</return>
            virtual std::size_t get_num_components_dynamic() const noexcept = 0;

            /// <summary>
            /// Get data array
            /// </summary>
            /// <return>Data array</return>
            virtual const std::vector<double>& get_data_dynamic() const = 0;

        protected:
            /// <summary>
            /// Default destructor
            /// </summary>
            virtual ~array_base() = default;

            /// Array name
            std::string name;

            /// Number of stored elements
            std::size_t num_elements;
        };
    }
}

#include "tpf_array_base.inl"
