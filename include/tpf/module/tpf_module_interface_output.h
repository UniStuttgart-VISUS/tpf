#pragma once

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Input interface for modules
        /// </summary>
        template <typename... Outputs>
        class interface_output
        {
        public:
            /// State that set_output(...) needs to be called before running the module
            static constexpr bool needs_output() { return true; }

            /// <summary>
            /// Constructor.
            /// </summary>
            interface_output();

            /// <summary>
            /// Set the output needed for algorithm execution.
            /// </summary>
            void set_output(Outputs... outputs);

        protected:
            /// <summary>
            /// Set the output needed for algorithm execution.
            /// The output types are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_output(Outputs... outputs) = 0;

            /// <summary>
            /// Perform a sanity check on the output
            /// </summary>
            void sanity_check() const;

        private:
            /// Information for sanity checks
            bool output_set;
        };
    }
}

#include "tpf_module_interface_output.inl"
