#pragma once

#include <exception>

namespace tpf
{
    namespace exception
    {
        /// <summary>
        /// Exception for handling trial-and-error in algorithms
        /// </summary>
        class trial_and_error_exception : public std::exception
        {
        public:
            trial_and_error_exception();
        };
    }
}

#include "tpf_trial_and_error_exception.inl"