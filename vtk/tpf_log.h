#pragma once

#include "../log/tpf_log.h"
#include "../log/tpf_log_base.h"

#include <functional>
#include <memory>
#include <string>

/// Define for creating a logger using ParaView's window or outstream objects (default: cerr and clog) for logging
#define __tpf_create_vtk_log() std::shared_ptr<tpf::log::log_base> tpf::log::logger = tpf::vtk::create_log_instance();

namespace tpf
{
    namespace vtk
    {
        /// <summary>
        /// Class for encapsulating VTK logging functions
        /// </summary>
        class vtk_log
        {
        public:
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="func">Function to call</param>
            /// <param name="active">Active stream?</param>
            vtk_log(std::function<void(const char*)> func, bool active = true);

            /// <summary>
            /// Stream operator
            /// </summary>
            friend void operator<<(vtk_log& stream, const std::string& input);

            /// <summary>
            /// Dummy flush
            /// </summary>
            void flush() const;

        private:
            /// Function
            std::function<void(const char*)> func;

            /// Active?
            bool active;
        };

        /// <summary>
        /// Create singleton instance for logging in VTK
        /// </summary>
        static std::shared_ptr<log::log_base> create_log_instance();
    }
}

#include "tpf_log.inl"