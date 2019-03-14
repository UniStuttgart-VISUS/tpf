#include "tpf_array_base.h"

#include <string>
#include <utility>

namespace tpf
{
    namespace data
    {
        inline array_base::array_base(std::string name, const std::size_t num_elements) noexcept
        {
            std::swap(this->name, name);
            this->num_elements = num_elements;
        }

        inline const std::string& array_base::get_name() const noexcept
        {
            return this->name;
        }

        inline void array_base::set_name(std::string name) noexcept
        {
            std::swap(this->name, name);
        }

        inline std::size_t array_base::get_num_elements() const noexcept
        {
            return this->num_elements;
        }
    }
}