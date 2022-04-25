#pragma once

#include "tpf_particle_seed.h"
#include "tpf_stream_trace.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/math/tpf_quaternion.h"

#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            /// <summary>
            /// Class for storing particle traces, original particle traces, rotations and seeds
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class path_trace : public stream_trace<floatp_t>
            {
            public:
                /// <summary>
                /// Disallow use of copy constructor
                /// </summary>
                path_trace(const path_trace&) = delete;

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source object</param>
                path_trace(path_trace&& move);

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source particle seed object to create this from</param>
                /// <param name="property_names">Name of attached property arrays</param>
                /// <param name="properties">Properties interpolated at the respective seed</param>
                path_trace(particle_seed<floatp_t>&& move, std::vector<std::string>&& property_names,
                    std::vector<std::vector<std::vector<double>>>&& properties);

                /// <summary>
                /// Return last original particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last original particle point</returns>
                const tpf::geometry::point<floatp_t>& get_last_original_particle(std::size_t index) const;

                /// <summary>
                /// Return last original particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last original particle point</returns>
                tpf::geometry::point<floatp_t>& get_last_original_particle(std::size_t index);

                /// <summary>
                /// Return last rotation for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last rotation</returns>
                const tpf::math::quaternion<floatp_t>& get_last_rotation(std::size_t index) const;

                /// <summary>
                /// Return last rotation for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last rotation</returns>
                tpf::math::quaternion<floatp_t>& get_last_rotation(std::size_t index);

                /// <summary>
                /// Add particle, original particle and rotation to a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle">New particle</param>
                /// <param name="original_particle">New original particle</param>
                /// <param name="sampled_velocity">Sampled velocity used to calculate the new particle's position</param>
                /// <param name="original_sampled_velocity">Sampled velocity used to calculate the new original particle's position</param>
                /// <param name="rotation">New rotation</param>
                /// <param name="properties">Properties interpolated at the particle</param>
                void add_particle(std::size_t index, const tpf::geometry::point<floatp_t>& particle,
                    const tpf::geometry::point<floatp_t>& original_particle, const Eigen::Matrix<floatp_t, 3, 1>& sampled_velocity,
                    const Eigen::Matrix<floatp_t, 3, 1>& original_sampled_velocity,
                    const tpf::math::quaternion<floatp_t>& rotation, std::vector<std::vector<double>>&& properties);

                /// <summary>
                /// Sort the particle traces for length and update the number of valid particles
                /// </summary>
                /// <param name="num_valid_particles">Previous number of valid particles</param>
                /// <returns>New number of valid particles</returns>
                virtual std::size_t sort_and_count(std::size_t num_valid_particles);

                /// <summary>
                /// Return the original sampled velocities for all particles
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Original sampled velocities</returns>
                const std::vector<Eigen::Matrix<floatp_t, 3, 1>>& get_original_sampled_velocities(std::size_t index) const;

            protected:
                /// Original particles
                std::vector<std::vector<tpf::geometry::point<floatp_t>>> original_particles;
                std::vector<std::vector<Eigen::Matrix<floatp_t, 3, 1>>> original_sampled_velocities;

                /// Rotations
                std::vector<std::vector<tpf::math::quaternion<floatp_t>>> rotations;
            };
        }
    }
}

#include "tpf_path_trace.inl"
