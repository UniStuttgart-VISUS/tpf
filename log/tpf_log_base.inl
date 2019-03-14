#include "tpf_log_base.h"

#include <mutex>
#include <string>

namespace tpf
{
    namespace log
    {
        inline void log_base::error_message(const std::string& message, const bool root_only) const
        {
            std::lock_guard<std::mutex> func_lock(this->lock);
            
            write_error_message(message, root_only);
        }

        inline void log_base::warning_message(const std::string& message, const bool root_only) const
        {
            std::lock_guard<std::mutex> func_lock(this->lock);
            
            write_warning_message(message, root_only);
        }

        inline void log_base::info_message(const std::string& message, const bool root_only) const
        {
            std::lock_guard<std::mutex> func_lock(this->lock);
            
            write_info_message(message, root_only);
        }
    }
}