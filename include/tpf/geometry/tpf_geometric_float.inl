#include "tpf_geometric_float.h"

#include <CGAL/number_utils.h>

#include <complex>
#include <type_traits>

namespace tpf
{
    namespace geometry
    {
        template <typename floatp_t, typename exact_t>
        inline geometric_float<floatp_t, exact_t>::geometric_float(exact_t value) : value(std::move(value))
        { }

        template <typename floatp_t, typename exact_t>
        inline geometric_float<floatp_t, exact_t>& geometric_float<floatp_t, exact_t>::operator=(exact_t value)
        {
            std::swap(this->value, std::move(value));

            return *this;
        }

        template <typename floatp_t, typename exact_t>
        inline geometric_float<floatp_t, exact_t>::operator exact_t() const
        {
            return this->value;
        }

        template <typename floatp_t, typename exact_t>
        inline exact_t geometric_float<floatp_t, exact_t>::get_exact_value() const
        {
            return this->value;
        }

        template <typename floatp_t, typename exact_t>
        template <typename F, typename E>
        inline geometric_float<floatp_t, exact_t>::operator typename std::enable_if<!std::is_same<F, E>::value, F>::type () const
        {
            return static_cast<floatp_t>(CGAL::to_double(this->value));
        }

        template <typename floatp_t, typename exact_t>
        inline floatp_t geometric_float<floatp_t, exact_t>::get_float_value() const
        {
            return static_cast<floatp_t>(CGAL::to_double(this->value));
        }

        template <typename floatp_t, typename exact_t>
        inline geometric_float<floatp_t, exact_t>::operator std::complex<floatp_t>() const
        {
            return static_cast<floatp_t>(CGAL::to_double(this->value));
        }

        template <typename floatp_t, typename exact_t>
        inline std::complex<floatp_t> geometric_float<floatp_t, exact_t>::get_complex_float_value() const
        {
            return static_cast<floatp_t>(CGAL::to_double(this->value));
        }
    }
}