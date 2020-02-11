#pragma once

#include <chrono>
#include <string>

namespace tpf
{
    namespace performance
    {
        /// <summary>
        /// Duration units
        /// </summary>
        /// <template name="duration_t">Chronos duration specialization</param>
        template <typename duration_t> struct duration_unit {};
        
        template <> struct duration_unit<std::chrono::nanoseconds> { static std::string value() { return std::string("ns"); } };
        template <> struct duration_unit<std::chrono::microseconds> { static std::string value() { return std::string("us"); } };
        template <> struct duration_unit<std::chrono::milliseconds> { static std::string value() { return std::string("ms"); } };
        template <> struct duration_unit<std::chrono::seconds> { static std::string value() { return std::string("s"); } };
        template <> struct duration_unit<std::chrono::minutes> { static std::string value() { return std::string("min"); } };
        template <> struct duration_unit<std::chrono::hours> { static std::string value() { return std::string("h"); } };
    }
}