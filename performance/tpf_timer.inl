#include "tpf_timer.h"

#include "tpf_duration_unit.h"

#include "../log/tpf_log.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace tpf
{
    namespace performance
    {
        template <typename duration_t>
        inline timer<duration_t>::timer(std::string name) : name(std::string(" (") + name + std::string(")")), start(std::chrono::steady_clock::now()), stopped(false)
        {
            log::info_message(__tpf_info_message("Timer started at ", timestamp(), this->name));
        }

        template <typename duration_t>
        inline timer<duration_t>::~timer()
        {
            if (!this->stopped)
            {
                stop();
            }
        }

        template <typename duration_t>
        inline void timer<duration_t>::stop()
        {
            this->stopped = true;

            const auto duration = std::chrono::duration_cast<duration_t>(std::chrono::steady_clock::now() - this->start);
            
            log::info_message(__tpf_info_message("Timer stopped", this->name, ": ", duration.count(), duration_unit<duration_t>::value()));
        }

        template <typename duration_t>
        inline std::string timer<duration_t>::timestamp() const
        {
            std::time_t now = std::time(nullptr);
            
            std::stringstream ss;
            ss << std::put_time(std::localtime(&now), "%X");
            return ss.str();
        }
    }
}