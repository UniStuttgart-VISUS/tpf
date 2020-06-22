#pragma once

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Parameter interface for modules
        /// </summary>
        template <typename... Parameters>
        class interface_parameters
        {
        public:
            /// State that set_parameters(...) needs to be called before running the module
            static constexpr bool needs_parameters() { return true; }

            /// <summary>
            /// Constructor.
            /// </summary>
            interface_parameters();

            /// <summary>
            /// Set the parameters needed for algorithm execution.
            /// </summary>
            void set_parameters(Parameters... parameters);

        protected:
            /// <summary>
            /// Set the parameters needed for algorithm execution.
            /// </summary>
            virtual void set_algorithm_parameters(Parameters... parameters) = 0;

            /// <summary>
            /// Perform a sanity check on the parameters
            /// </summary>
            void sanity_check() const;

        private:
            /// Information for sanity checks
            bool parameters_set;
        };
    }
}

#include "tpf_module_interface_parameters.inl"
