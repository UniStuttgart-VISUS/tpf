#pragma once

#include <complex>
#include <type_traits>

namespace tpf
{
    namespace geometry
    {
        /// <summary>
        /// Stores an exact value and allows for implicit conversion into a floating point value
        /// </summary>
        /// <template name="floatp_t">Floating point type</template>
        /// <template name="exact_t">Exact value type</template>
        template <typename floatp_t, typename exact_t>
        class geometric_float
        {
        public:
            /// Value and floating point types
            using value_type = exact_t;
            using float_type = floatp_t;

            /// <summary>
            /// Create or assign from an exact value
            /// </summary>
            /// <param name="value">Exact value</param>
            geometric_float(exact_t value);
            geometric_float& operator=(exact_t value);

            /// <summary>
            /// Default copy/move constructor/assignment
            /// </summary>
            geometric_float(const geometric_float&) = default;
            geometric_float(geometric_float&&) noexcept = default;
            geometric_float& operator=(const geometric_float&) = default;
            geometric_float& operator=(geometric_float&&) noexcept = default;

            /// <summary>
            /// Return exact value
            /// </summary>
            /// <returns>Stored exact value</returns>
            operator exact_t() const;
            exact_t get_exact_value() const;

            /// <summary>
            /// Return inexact floating point value
            /// </summary>
            /// <returns>Inexact floating point value</returns>
            operator floatp_t() const;
            floatp_t get_float_value() const;

            operator std::complex<floatp_t>() const;
            std::complex<floatp_t> get_complex_float_value() const;

        private:
            /// Stored value
            exact_t value;
        };
    }
}

#include "tpf_geometric_float.inl"
