#pragma once

#include "tpf_generic_module_base.h"

#include <functional>
#include <optional>
#include <string>
#include <type_traits>

namespace tpf
{
    namespace modules
    {
        /// <summary>
        /// Module base class
        /// </summary>
        template <typename Input, typename Output, typename Parameters, typename Callbacks>
        class module_base : public generic_module_base
        {
        public:
            using input_t = Input;
            using output_t = Output;
            using parameters_t = Parameters;
            using callbacks_t = Callbacks;

            static constexpr bool needs_input() { return !std::is_same<input_t, void>::value; }
            static constexpr bool needs_output() { return !std::is_same<output_t, void>::value; }
            static constexpr bool needs_parameters() { return !std::is_same<parameters_t, void>::value; }
            static constexpr bool needs_callbacks() { return !std::is_same<callbacks_t, void>::value; }

        public:
            /// <summary>
            /// Constructor.
            /// Derived class has to state need for setting input, output, parameters, and callbacks.
            /// </summary>
            module_base();

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run() override final;

            /// <summary>
            /// Set the input needed for algorithm execution.
            /// The input types are defined by the derived class.
            /// </summary>
            template <typename in_t = Input>
            void set_input(typename std::enable_if<!std::is_same<in_t, void>::value, input_t>::type inputs);

            template <typename in_t = Input>
            void set_input(typename std::enable_if<std::is_same<in_t, void>::value>::type);

            /// <summary>
            /// Set the output resulting from algorithm execution.
            /// The output types are defined by the derived class.
            /// </summary>
            template <typename out_t = Output>
            void set_output(typename std::enable_if<!std::is_same<out_t, void>::value, output_t>::type outputs);

            template <typename out_t = Output>
            void set_output(typename std::enable_if<std::is_same<out_t, void>::value>::type);

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
            /// Set the input needed for algorithm execution.
            /// The input types are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_input(input_t) = 0;

            /// <summary>
            /// Set the output resulting from algorithm execution.
            /// The output types are defined by the derived class.
            /// </summary>
            virtual void set_algorithm_output(output_t) = 0;

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

            /// <summary>
            /// Process an optional parameter and set its default value if none was given
            /// </summary>
            template <typename T>
            T get_or_default(std::optional<T> parameter, T default_value) const;

            /// <summary>
            /// Process an optional reference parameter and set its default value if none was given
            /// </summary>
            template <typename T>
            T* get_or_default(std::optional<std::reference_wrapper<T>> parameter, T* default_value = nullptr) const;

        private:
            /// Information for sanity checks
            bool input_set;
            bool output_set;
            bool parameters_set;
            bool callbacks_set;
        };
    }
}

#include "tpf_module_base.inl"
