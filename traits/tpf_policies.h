#pragma once

#include "tpf_stdext.h"

#include "../policies/tpf_interpolatable.h"

namespace tpf
{
    namespace traits
    {        
        /// <summary>
        /// Checks if a data structure supports interpolation
        /// </summary>
        /// <template name=""></template>
        template <typename class_t>
        struct supports_interpolatable
        {
            static constexpr bool value = is_base_of_template<policies::interpolatable, class_t>::value && class_t::value_type_interpolatable;
        };
    }
}