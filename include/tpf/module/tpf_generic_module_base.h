#pragma once

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Generic module base class to avoid template parameters in the base class
        /// </summary>
        class generic_module_base
        {
        public:
            /// <summary>
            /// Run module
            /// </summary>
            virtual void run() = 0;
        };
    }
}
