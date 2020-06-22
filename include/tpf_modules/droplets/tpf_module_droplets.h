#pragma once

#include "tpf/module/tpf_module_base.h"
#include "tpf/module/tpf_module_interface_input.h"
#include "tpf/module/tpf_module_interface_output.h"
#include "tpf/module/tpf_module_interface_parameters.h"

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
            interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>,
            interface_output<data::grid<long long, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
                std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>, data::polydata<float_t>&>,
            interface_parameters<bool, bool, bool, bool, droplets_aux::scale_method_t, droplets_aux::rotation_method_t>>
        {
            using input_t = interface_input<const data::grid<float_t, float_t, 3, 1>&, const data::grid<float_t, float_t, 3, 3>&,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>>>;
            using output_t = interface_output<data::grid<long long, float_t, 3, 1>&, data::grid<float_t, float_t, 3, 1>&,
                std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>>, data::polydata<float_t>&>;
            using parameters_t = interface_parameters<bool, bool, bool, bool, droplets_aux::scale_method_t, droplets_aux::rotation_method_t>;

            using base_t = module_base<input_t, output_t, parameters_t>;

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
            /// <param name="fractions">Input fractions</param>
            /// <param name="positions">Input positions</param>
            /// <param name="velocities">Input velocities</param>
            virtual void set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions, const data::grid<float_t, float_t, 3, 3>& positions,
                std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>> velocities) override;

            /// <summary>
            /// Set output
            /// </summary>
            /// <param name="droplet_ids">Output droplet IDs</param>
            /// <param name="droplet_volumes">Output volumes</param>
            /// <param name="droplet_velocities">Output velocities</param>
            /// <param name="droplets">Output droplets</param>
            virtual void set_algorithm_output(data::grid<long long, float_t, 3, 1>& droplet_ids, data::grid<float_t, float_t, 3, 1>& droplet_volumes,
                std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>> droplet_velocities, data::polydata<float_t>& droplets) override;

            /// <summary>
            /// Set parameter
            /// </summary>
            /// <param name="calculate_translation">Calculate translational velocity</param>
            /// <param name="calculate_rotation">Calculate rotational velocity</param>
            /// <param name="calculate_energy">Calculate energy</param>
            /// <param name="calculate_inertia">Calculate inertia</param>
            /// <param name="scaling_method">Scaling method</param>
            /// <param name="rotation_method">Method for computing rotation</param>
            virtual void set_algorithm_parameters(bool calculate_translation, bool calculate_rotation, bool calculate_energy, bool calculate_inertia,
                droplets_aux::scale_method_t scaling_method, droplets_aux::rotation_method_t rotation_method) override;

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
