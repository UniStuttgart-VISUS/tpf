#pragma once

#include "tpf/module/tpf_module_base.h"

#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include <functional>
#include <optional>
#include <string>
#include <tuple>

namespace tpf
{
    namespace modules
    {
        namespace droplets_aux
        {
            /// <summary>
            /// Scaling method
            /// </summary>
            enum class scale_method_t
            {
                INVALID = -1,
                NORMALIZED, BY_EIGENVALUE, BY_LOGARITHMIC_EIGENVALUE
            };

            /// <summary>
            /// Rotation method
            /// </summary>
            enum class rotation_method_t
            {
                INVALID = -1,
                MECHANICS, VELOCITIES, PCA
            };
        }

        /// <summary>
        /// Module for calculation of droplet IDs
        /// </summary>
        /// <template name="float_t">Floating point type</template>
        template <typename float_t>
        class droplets : public module_base<
            // Input
            std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
            std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>,
            // Output
            std::tuple<data::grid<long long, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
            std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>, data::polydata<float_t>&>,
            // Parameters
            std::tuple<bool, bool, bool, bool, droplets_aux::scale_method_t, droplets_aux::rotation_method_t>,
            // Callbacks
            void>
        {
            using base_t = module_base<std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>, std::tuple<data::grid<long long, float_t, 3, 1>&,
                data::grid<float_t, float_t, 3, 1>&, std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>, data::polydata<float_t>&>,
                std::tuple<bool, bool, bool, bool, droplets_aux::scale_method_t, droplets_aux::rotation_method_t>, void>;

        public:
            /// <summary>
            /// Return the number of ghost levels
            /// </summary>
            /// <returns>Number of ghost levels</returns>
            static std::size_t get_num_required_ghost_levels();

            /// <summary>
            /// Constructor
            /// </summary>
            droplets();

            /// <summary>
            /// Return the name of the module
            /// </summary>
            /// <returns>Module name</returns>
            std::string get_name() const override;

        protected:
            /// <summary>
            /// Set input
            /// </summary>
            /// <param name="input">Input [fractions, positions, velocities]</param>
            virtual void set_algorithm_input(std::tuple<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>> input) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="output">Output [droplet IDs, volumes, velocities, droplets]</param>
            virtual void set_algorithm_output(std::tuple<data::grid<long long, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
                std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>, data::polydata<float_t>&> output) override;

            /// <summary>
            /// Set parameter
            /// </summary>
            /// <param name="parameters">[Calculate translational velocity, Calculate rotational velocity,
            /// Calculate energy, Calculate inertia, Scaling method, Method for computing rotation]</param>
            virtual void set_algorithm_parameters(std::tuple<bool, bool, bool, bool,
                droplets_aux::scale_method_t, droplets_aux::rotation_method_t> parameters) override;

            /// <summary>
            /// Run module
            /// </summary>
            void run_algorithm() override;

        private:
            /// <summary>
            /// Scale value
            /// </summary>
            /// <param name="value">Value</param>
            /// <returns>Scaled value</returns>
            float_t scale(float_t value) const;

            /// Fractions
            const data::grid<float_t, float_t, 3, 1>* fractions;

            /// Velocities
            const data::grid<float_t, float_t, 3, 3>* velocities;

            /// Positions
            const data::grid<float_t, float_t, 3, 3>* positions;

            /// Droplet IDs
            data::grid<long long, float_t, 3, 1>* droplet_ids;

            /// Volumes
            data::grid<float_t, float_t, 3, 1>* droplet_volumes;

            /// Global velocities
            data::grid<float_t, float_t, 3, 3>* droplet_velocities;

            /// Droplets
            data::polydata<float_t>* local_droplets;

            /// Deciders on what to calculate
            bool calculate_translation;
            bool calculate_rotation;
            bool calculate_energy;
            bool calculate_inertia;

            /// Scaling method
            droplets_aux::scale_method_t scale_method;

            /// Method for computing rotation
            droplets_aux::rotation_method_t rotation_method;
        };
    }
}

#include "tpf_module_droplets.inl"
