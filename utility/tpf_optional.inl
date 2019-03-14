#include "tpf_optional.h"

#include <exception>
#include <memory>
#include <utility>

namespace tpf
{
    namespace utility
    {
        template <typename t>
        inline optional<t>::optional() : valid(false), value(nullptr)
        { }

        template <typename t>
        inline optional<t>::optional(nullopt_t) : valid(false), value(nullptr)
        { }

        template <typename t>
        template <class u>
        inline optional<t>::optional(u&& value) : valid(true)
        {
            this->value = std::make_shared<t>(value);
        }

        template <typename t>
        inline optional<t>& optional<t>::operator=(nullopt_t)
        {
            this->value = nullptr;
            this->valid = false;
        }

        template <typename t>
        template <class u>
        inline optional<t>& optional<t>::operator=(u&& value)
        {
            this->value = std::make_shared<t>(value);
            this->valid = true;
        }

        template <typename t>
        inline std::shared_ptr<const t> optional<t>::operator->() const
        {
            if (this->valid)
            {
                return this->value;
            }
            else
            {
                return nullptr;
            }
        }

        template <typename t>
        inline std::shared_ptr<t> optional<t>::operator->()
        {
            if (this->valid)
            {
                return this->value;
            }
            else
            {
                return nullptr;
            }
        }

        template <typename t>
        inline const t& optional<t>::operator*() const&
        {
            if (this->valid)
            {
                return *this->value;
            }
            else
            {
                throw std::exception();
            }
        }

        template <typename t>
        inline t& optional<t>::operator*() &
        {
            if (this->valid)
            {
                return *this->value;
            }
            else
            {
                throw std::exception();
            }
        }

        template <typename t>
        inline const t&& optional<t>::operator*() const&&
        {
            if (this->valid)
            {
                return std::move(*this->value);
            }
            else
            {
                throw std::exception();
            }
        }

        template <typename t>
        inline t&& optional<t>::operator*() &&
        {
            if (this->valid)
            {
                return std::move(*this->value);
            }
            else
            {
                throw std::exception();
            }
        }

        template <typename t>
        inline optional<t>::operator bool() const
        {
            return this->valid;
        }
    }
}