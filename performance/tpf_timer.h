#pragma once

#include <chrono>
#include <string>

#ifdef __tpf_performance
#define __tpf_timer(name, duration) tpf::performance::timer<duration> __tpf_my_timer__(name);
#define __tpf_timer_start(var, name, duration) tpf::performance::timer<duration> var(name);
#define __tpf_timer_stop(var) var.stop();
#else
#define __tpf_timer(name, duration)
#define __tpf_timer_start(var, name, duration)
#define __tpf_timer_stop(var)
#endif

namespace tpf
{
    namespace performance
    {
        /// <summary>
        /// Timer class for performance measures
        /// </summary>
        /// <template name="duration_t">Duration</template>
        template <typename duration_t = std::chrono::milliseconds>
        class timer
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="name">Timer name for logging</param>
            timer(std::string name = std::string(""));

            /// <summary>
            /// Destructor stops the timer
            /// </summary>
            virtual ~timer();
            
            /// <summary>
            /// Stop the timer and output measured duration
            /// </summary>
            void stop();
            
        private:
            /// <summary>
            /// Create timestamp of the current system time
            /// </summary>
            /// <returns>Timestamp of now</returns>
            std::string timestamp() const;

            /// Timer name
            const std::string name;
            
            /// Start of the timer
            const std::chrono::time_point<std::chrono::steady_clock> start;
            
            /// Has the timer already been stopped manually?
            bool stopped;
        };
    }
}

#include "tpf_timer.inl"