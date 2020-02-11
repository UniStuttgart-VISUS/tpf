#include "tpf_log.h"

#include "tpf_log_base.h"
#include "tpf_log_ostream.h"

#include <memory>
#include <sstream>
#include <string>

namespace tpf
{
    namespace log
    {
        template <typename Err_Stream, typename Warn_Stream, typename Info_Stream>
        inline std::shared_ptr<log_base> create_instance(Err_Stream& err_stream, Warn_Stream& warn_stream, Info_Stream& info_stream)
        {
            return std::dynamic_pointer_cast<log_base>(std::make_shared<log_ostream<Err_Stream, Warn_Stream, Info_Stream>>(err_stream, warn_stream, info_stream));
        }

        inline void error_message(const std::string& message, const bool root_only)
        {
            return logger->error_message(message, root_only);
        }

        inline void warning_message(const std::string& message, const bool root_only)
        {
            return logger->warning_message(message, root_only);
        }

        inline void info_message(const std::string& message, const bool root_only)
        {
            return logger->info_message(message, root_only);
        }

        namespace
        {
            template <typename First>
            inline std::string concatenate_message(const First& first)
            {
                std::stringstream ss;
                ss << first;

                return ss.str();
            }

            template <typename First, typename... Rest>
            inline std::string concatenate_message(const First& first, const Rest&... rest)
            {
                std::stringstream ss;
                ss << first << concatenate_message(rest...);

                return ss.str();
            }
        }

        template <typename... MessageT>
        inline std::string create_message(const MessageT&... message)
        {
            return concatenate_message(message...);
        }
    }
}