#pragma once

#include <stdexcept>

namespace tpf
{
    namespace exception
    {
        /// <summary>
        /// Exception for handling functions which are not yet implemented
        /// </summary>
        class not_implemented_exception : public std::logic_error
        {
        public:
            explicit not_implemented_exception(std::string message = "");
        };
    }
}

#include "tpf_not_implemented_exception.inl"