#pragma once

#include "tpf_log_base.h"

#include <iostream>
#include <ostream>
#include <memory>
#include <string>

/// Define for creating a logger using outstream objects (default: cerr and clog) for logging
#define __tpf_create_log(...) std::shared_ptr<tpf::log::log_base> tpf::log::logger = tpf::log::create_instance(__VA_ARGS__);

/// Defines for creating error, warning and information messages
#ifdef __tpf_log_prefix
#define __tpf_error_message(...) ::tpf::log::create_message("[ERROR]: ", __VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#define __tpf_warning_message(...) ::tpf::log::create_message("[WARNING]: ", __VA_ARGS__)
#define __tpf_info_message(...) ::tpf::log::create_message("[INFO]: ", __VA_ARGS__)

#ifdef __tpf_debug
#undef __tpf_warning_message
#undef __tpf_info_message

#define __tpf_warning_message(...) ::tpf::log::create_message("[WARNING]: ", __VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#define __tpf_info_message(...) ::tpf::log::create_message("[INFO]: ", __VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#endif
#else
#define __tpf_error_message(...) ::tpf::log::create_message(__VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#define __tpf_warning_message(...) ::tpf::log::create_message(__VA_ARGS__)
#define __tpf_info_message(...) ::tpf::log::create_message(__VA_ARGS__)

#ifdef __tpf_debug
#undef __tpf_warning_message
#undef __tpf_info_message

#define __tpf_warning_message(...) ::tpf::log::create_message(__VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#define __tpf_info_message(...) ::tpf::log::create_message(__VA_ARGS__, "\n  File: ", __FILE__, ", line: ", __LINE__)
#endif
#endif

/// Defines for creating nested messages
#define __tpf_nested_error_message(inner, ...) __tpf_error_message(__VA_ARGS__, "\n ", inner)
#define __tpf_nested_warning_message(inner, ...) __tpf_warning_message(__VA_ARGS__, "\n ", inner)
#define __tpf_nested_info_message(inner, ...) __tpf_info_message(__VA_ARGS__, "\n ", inner)

/// Defines for development debugging
#define __tpf_dev_debug_output tpf::log::info_message(__tpf_info_message(__FILE__, " ~ ", __LINE__));

namespace tpf
{
    namespace log
    {
        /// <summary>
        /// Create the singleton instance
        /// </summary>
        /// <param name="err_stream">Output error stream</param>
        /// <param name="warn_stream">Output warning stream</param>
        /// <param name="info_stream">Output info stream</param>
        /// <returns>Created instance</returns>
        template <typename Err_Stream = std::ostream, typename Warn_Stream = std::ostream, typename Info_Stream = std::ostream>
        static std::shared_ptr<log_base> create_instance(Err_Stream& err_stream = std::cerr, Warn_Stream& warn_stream = std::clog, Info_Stream& info_stream = std::clog);

        /// <summary>
        /// Write error message
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
        static void error_message(const std::string& message, bool root_only = false);

        /// <summary>
        /// Write warning message
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
        static void warning_message(const std::string& message, bool root_only = false);

        /// <summary>
        /// Write info message
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
        static void info_message(const std::string& message, bool root_only = false);

        /// <summary>
        /// Create a message
        /// </summary>
        /// <template name="MessageT">Type of the message</template>
        /// <param name="message">The message</param>
        /// <returns>Message</returns>
        template <typename... MessageT>
        static std::string create_message(const MessageT&... message);

        /// Instance for logging
        extern std::shared_ptr<log_base> logger;
    }
}

#include "tpf_log.inl"