#include "tpf_bool.h"

namespace tpf
{
    inline bool_t::bool_t()
    {}

    inline bool_t::bool_t(const bool value)
    {
        this->value = value;
    }

    inline bool_t::bool_t(const bool_t& value)
    {
        this->value = value.value;
    }

    inline bool_t::bool_t(bool_t&& value) noexcept
    {
        this->value = value.value;
    }

    inline bool_t& bool_t::operator=(const bool_t& value)
    {
        this->value = value.value;

        return *this;
    }

    inline bool_t& bool_t::operator=(bool_t&& value) noexcept
    {
        this->value = value.value;

        return *this;
    }

    inline bool_t& bool_t::operator=(bool value)
    {
        this->value = value;

        return *this;
    }

    inline bool_t::operator bool&()
    {
        return this->value;
    }

    inline bool_t::operator bool() const
    {
        return this->value;
    }
}
