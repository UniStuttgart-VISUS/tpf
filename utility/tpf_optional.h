#pragma once

#include <memory>

namespace tpf
{
    namespace utility
    {
        /// Struct indicating optional class to set validity to false
        struct nullopt_t {};
        static nullopt_t nullopt;

        /// <summary>
        /// Optionally store a value of given type
        /// </summary>
        /// <template name="t">Stored type</template>
        template <typename t>
        class optional
        {
        public:
            using value_type = t;

            /// <summary>
            /// Constructor
            /// </summary>
            optional();
            optional(nullopt_t);
            optional(const optional&) = default;
            optional(optional&&) = default;
            
            /// <summary>
            /// Store value and set to valid
            /// </summary>
            template <class u = t>
            optional(u&& value);

            /// <summary>
            /// Assign from other optional
            /// </summary>
            optional& operator=(nullopt_t);
            optional& operator=(const optional&) = default;
            optional& operator=(optional&&) = default;

            /// <summary>
            /// Store value and set to valid
            /// </summary>
            template <class u = t>
            optional& operator=(u&& value);

            /// <summary>
            /// Access stored value
            /// </summary>
            /// <returns>Stored value</returns>
            std::shared_ptr<const t> operator->() const;
            std::shared_ptr<t> operator->();
            const t& operator*() const&;
            t& operator*() &;
            const t&& operator*() const&&;
            t&& operator*() &&;

            /// <summary>
            /// Answers validity of stored value
            /// </summary>
            /// <returns>Validity</returns>
            explicit operator bool() const;

        private:
            /// Stored value
            std::shared_ptr<t> value;

            /// Is stored value valid?
            bool valid;
        };
    }
}

#include "tpf_optional.inl"