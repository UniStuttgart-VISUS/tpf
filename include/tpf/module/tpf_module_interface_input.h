#pragma once

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Input interface for modules
        /// </summary>
        template <typename... Inputs>
        class interface_input
        {
        public:
            /// State that set_input(...) needs to be called before running the module
            static constexpr bool needs_input() { return true; }

            /// <summary>
            /// Constructor.
            /// </summary>
            interface_input();

            /// <summary>
            /// Set the input needed for algorithm execution.
            /// </summary>
            void set_input(Inputs... inputs);

        protected:
            /// <summary>
            /// Set the input needed for algorithm execution.
            /// </summary>
            virtual void set_algorithm_input(Inputs... inputs) = 0;

            /// <summary>
            /// Perform a sanity check on the input
            /// </summary>
            void sanity_check() const;

        private:
            /// Information for sanity checks
            bool input_set;
        };
    }
}

#include "tpf_module_interface_input.inl"
