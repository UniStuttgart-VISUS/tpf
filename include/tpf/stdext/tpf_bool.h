#pragma once

namespace tpf
{
    /// <summary>
    /// Boolean type for use in, e.g., vectors
    /// </summary>
    struct bool_t
    {
    public:
        /// <summary>
        /// Constructors
        /// </summary>
        bool_t();
        bool_t(bool value);
        bool_t(const bool_t&) = default;
        bool_t(bool_t&&) noexcept = default;

        /// <summary>
        /// Assignments
        /// </summary>
        bool_t& operator=(const bool_t&) = default;
        bool_t& operator=(bool_t&&) noexcept = default;
        bool_t& operator=(bool value);

        /// <summary>
        /// Conversions
        /// </summary>
        operator bool&();
        operator bool() const;

    private:
        /// Stored value
        bool value;
    };
}

#include "tpf_bool.inl"
