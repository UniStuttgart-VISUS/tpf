#pragma once

#include "tpf_geometric_object.h"

#include "../utility/tpf_optional.h"

#include <memory>
#include <type_traits>
#include <vector>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Safely create the specified geometry
        /// </summary>
        /// <template name="object_t">Geometry object type</template>
        /// <template name="arguments_t">Argument types for the respective constructor</template>
        /// <param name="arguments_t">Arguments for the respective constructor</param>
        template <typename object_t, typename... arguments_t>
        static utility::optional<object_t> make(arguments_t... args) noexcept;

        /// <summary>
        /// Safely create the specified geometry
        /// </summary>
        /// <template name="object_t">Geometry object type</template>
        /// <template name="arguments_t">Argument types for the respective constructor</template>
        /// <param name="arguments_t">Arguments for the respective constructor</param>
        template <typename object_t, typename... arguments_t>
        static typename std::enable_if<!std::is_floating_point<object_t>::value, std::shared_ptr<object_t>>::type make_shared(arguments_t... args) noexcept;
        
        /// <summary>
        /// Safely create the specified geometry
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">CGAL kernel</template>
        /// <template name="arguments_t">Argument types for the respective constructor</template>
        /// <param name="geometry">Geometry type</param>
        /// <param name="arguments_t">Arguments for the respective constructor</param>
        template <typename floatp_t, typename kernel_t = default_kernel_t, typename... arguments_t>
        static typename std::enable_if<std::is_floating_point<floatp_t>::value, std::shared_ptr<geometric_object<floatp_t>>>::type
            make_shared(geometry_t geometry, arguments_t... args) noexcept;

        /// <summary>
        /// Deserialize object from a binary representation, preceded by an indicator for the object type
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="kernel_t">CGAL kernel</template>
        /// <param name="serialized">Serialized object</param>
        template <typename floatp_t, typename kernel_t = default_kernel_t>
        static std::shared_ptr<geometric_object<floatp_t>> deserialize(const std::vector<char>& serialized);
    }
}

#include "tpf_factory.inl"