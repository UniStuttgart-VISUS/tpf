#pragma once

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Input interface for modules
        /// </summary>
        template <typename... Callbacks>
        class interface_callbacks
        {
        public:
            /// State that set_callback(...) needs to be called before running the module
            static constexpr bool needs_callbacks() { return true; }

            /// <summary>
            /// Constructor.
            /// </summary>
            interface_callbacks();

            /// <summary>
            /// Set the callback needed for algorithm execution.
            /// </summary>
            void set_callbacks(Callbacks... callbacks);

        protected:
            /// <summary>
            /// Set the callback needed for algorithm execution.
            /// The callback types are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_callbacks(Callbacks... callbacks) = 0;

            /// <summary>
            /// Perform a sanity check on the callbacks
            /// </summary>
            void sanity_check() const;

        private:
            /// Information for sanity checks
            bool callbacks_set;
        };
    }
}

#include "tpf_module_interface_callbacks.inl"
