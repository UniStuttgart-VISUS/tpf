#pragma once

#include "tpf_generic_module_base.h"

#include <functional>
#include <optional>
#include <string>

namespace tpf
{
    namespace modules
    {
        template <typename T> using opt_arg = std::optional<std::reference_wrapper<T>>;

        /// <summary>
        /// Module base class
        /// </summary>
        template <typename... Interfaces>
        class module_base : public generic_module_base, public Interfaces...
        {
        public:
            /// <summary>
            /// Constructor.
            /// </summary>
            module_base();

            /// <summary>
            /// Run module
            /// </summary>
            virtual void run() override final;

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
            /// Perform a sanity check on input, output, parameters and callbacks
            /// </summary>
            void sanity_checks() const;
        };

        /// <summary>
        /// Process an optional parameter and set its default value if none was given
        /// </summary>
        template <typename T>
        T get_or_default(std::optional<T> parameter, T default_value);

        /// <summary>
        /// Process an optional reference parameter and set its default value if none was given
        /// </summary>
        template <typename T>
        T* get_or_default(std::optional<std::reference_wrapper<T>> parameter, T* default_value = nullptr);

        /// <summary>
        /// Set an optional parameter or provide nothing
        /// </summary>
        template <typename T>
        std::optional<std::reference_wrapper<T>> set_or_default(T& parameter, bool set);

        /// <summary>
        /// Set an optional parameter or provide nothing
        /// </summary>
        template <typename T>
        std::optional<std::reference_wrapper<T>> set_or_default(T* parameter);
    }
}

#include "tpf_module_base.inl"
