#pragma once

#include <mutex>
#include <string>

namespace tpf
{
    namespace log
    {
        /// <summary>
        /// Base class for logging
        /// </summary>
        class log_base
        {
        public:
            /// <summary>
            /// Write error message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            void error_message(const std::string& message, bool root_only = false) const;

            /// <summary>
            /// Write warning message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            void warning_message(const std::string& message, bool root_only = false) const;

            /// <summary>
            /// Write info message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            void info_message(const std::string& message, bool root_only = false) const;

        private:
            /// <summary>
            /// Write error message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_error_message(const std::string& message, bool root_only) const = 0;

            /// <summary>
            /// Write warning message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_warning_message(const std::string& message, bool root_only) const = 0;

            /// <summary>
            /// Write info message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_info_message(const std::string& message, bool root_only) const = 0;

            /// Lock for proper output in asynchronous environments
            mutable std::mutex lock;
        };
    }
}

#include "tpf_log_base.inl"