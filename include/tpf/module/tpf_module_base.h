#pragma once

#include <limits>
#include <string>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module base class
        /// </summary>
        template <typename Parameters, typename Callbacks>
        class module_base
        {
        public:
            using parameters_t = Parameters;
            using callbacks_t = Callbacks;

            /// <summary>
            /// Constructor
            /// </summary>
            module_base();

            /// <summary>
            /// Run module
            /// </summary>
            void run();

            /// <summary>
            /// Set the parameters needed for algorithm execution.
            /// The parameter types are defined by the derived class.
            /// </summary>
            template <typename param_t = Parameters>
            void set_parameters(typename std::enable_if<!std::is_same<param_t, void>::value, parameters_t>::type parameters);

            template <typename param_t = Parameters>
            void set_parameters(typename std::enable_if<std::is_same<param_t, void>::value>::type);

            /// <summary>
            /// Set the callback functions needed for algorithm execution.
            /// The function signatures are defined by the derived class.
            /// </summary>
            template <typename cb_t = Callbacks>
            void set_callbacks(typename std::enable_if<!std::is_same<cb_t, void>::value, callbacks_t>::type callbacks);

            template <typename cb_t = Callbacks>
            void set_callbacks(typename std::enable_if<std::is_same<cb_t, void>::value>::type);

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            virtual std::string get_name() const = 0;

        protected:
            /// <summary>
            /// Run module's algorithm
            /// </summary>
            virtual void run_algorithm() = 0;

            /// <summary>
            /// Set the parameters needed for algorithm execution.
            /// The parameter types are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_parameters(parameters_t);

            /// <summary>
            /// Set the callback functions needed for algorithm execution.
            /// The function signatures are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_callbacks(callbacks_t);

            /// <summary>
            /// Perform a sanity check on set parameters and callbacks
            /// </summary>
            virtual bool sanity_check() const;

        private:
            /// Information for sanity checks
            bool parameters_set;
            bool callbacks_set;
        };
    }
}

#include "tpf_module_base.inl"