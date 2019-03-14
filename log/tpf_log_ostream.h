#pragma once

#include "tpf_log_base.h"

#include "../traits/tpf_function.h"

#ifdef __tpf_use_mpi
#include "mpi.h"
#endif

#include <ostream>
#include <string>

namespace tpf
{
    namespace log
    {
        /// <summary>
        /// Default implementation for logging using outstream objects
        /// </summary>
        /// <template name="error_stream_t">Error output stream type</template>
        /// <template name="warning_stream_t">Warning output stream type</template>
        /// <template name="info_stream_t">Info output stream type</template>
        template <typename error_stream_t = std::ostream, typename warning_stream_t = std::ostream, typename info_stream_t = std::ostream>
        class log_ostream : public log_base
        {
            __tpf_has_function(flush);

            static_assert(has_function_flush<error_stream_t>::value, "Error stream must implement function 'flush'.");
            static_assert(has_function_flush<warning_stream_t>::value, "Warning stream must implement function 'flush'.");
            static_assert(has_function_flush<info_stream_t>::value, "Info stream must implement function 'flush'.");

        public:
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="err_stream">Output error stream</param>
            /// <param name="warn_stream">Output warning stream</param>
            /// <param name="info_stream">Output info stream</param>
            log_ostream(error_stream_t& err_stream, warning_stream_t& warn_stream, info_stream_t& info_stream);

        protected:
            /// <summary>
            /// Write error message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_error_message(const std::string& message, bool root_only) const;

            /// <summary>
            /// Write warning message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_warning_message(const std::string& message, bool root_only) const;

            /// <summary>
            /// Write info message
            /// </summary>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            virtual void write_info_message(const std::string& message, bool root_only) const;

        private:
            /// <summary>
            /// 
            /// </summary>
            /// <template name="stream_t">Stream type</template>
            /// <param name="stream">Output stream</param>
            /// <param name="message">The message</param>
            /// <param name="root_only">Print the message only on the root node (MPI rank 0)</param>
            template <typename stream_t>
            void write_message(stream_t& stream, const std::string& message, bool root_only) const;

            /// Output streams
            error_stream_t& err_stream;
            warning_stream_t& warn_stream;
            info_stream_t& info_stream;

#ifdef __tpf_use_mpi
            /// MPI communicator
            mutable MPI_Comm world_comm;
#endif

            /// Process information
            mutable int rank;
            mutable int num_processes;
        };
    }
}

#include "tpf_log_ostream.inl"