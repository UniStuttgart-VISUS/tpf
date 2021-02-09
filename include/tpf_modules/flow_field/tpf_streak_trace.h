#pragma once

#include "tpf_particle_seed.h"
#include "tpf_path_trace.h"

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
            /// Class for storing particle traces, original particle traces, rotations, seeds and seed rotations
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class streak_trace : public path_trace<floatp_t>
            {
            public:
                /// <summary>
                /// Disallow use of copy constructor
                /// </summary>
                streak_trace(const streak_trace&) = delete;

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source object</param>
                streak_trace(streak_trace&& move);

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source particle seed object to create this from</param>
                streak_trace(particle_seed<floatp_t>&& move);

                /// <summary>
                /// Get number of particles in a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Number of particles in the trace</returns>
                std::size_t get_num_particles(std::size_t index) const;

                /// <summary>
                /// Return specific particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the particle in the trace</param>
                /// <returns>Particle point</returns>
                const tpf::geometry::point<floatp_t>& get_particle(std::size_t index, std::size_t particle_index) const;

                /// <summary>
                /// Return specific particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the particle in the trace</param>
                /// <returns>Particle point</returns>
                tpf::geometry::point<floatp_t>& get_particle(std::size_t index, std::size_t particle_index);

                /// <summary>
                /// Return specific original particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the original particle in the trace</param>
                /// <returns>Original particle point</returns>
                const tpf::geometry::point<floatp_t>& get_original_particle(std::size_t index, std::size_t particle_index) const;

                /// <summary>
                /// Return specific original particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the original particle in the trace</param>
                /// <returns>Original particle point</returns>
                tpf::geometry::point<floatp_t>& get_original_particle(std::size_t index, std::size_t particle_index);

                /// <summary>
                /// Return specific rotation for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the rotation in the trace</param>
                /// <returns>Rotation</returns>
                const tpf::math::quaternion<floatp_t>& get_rotation(std::size_t index, std::size_t particle_index) const;

                /// <summary>
                /// Return specific rotation for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle_index">Index of the rotation in the trace</param>
                /// <returns>Rotation</returns>
                tpf::math::quaternion<floatp_t>& get_rotation(std::size_t index, std::size_t particle_index);

                /// <summary>
                /// Return the specified original, advected seed point
                /// </summary>
                /// <param name="index">Index of the seed</param>
                /// <returns>Seed point</returns>
                const tpf::geometry::point<floatp_t>& get_original_seed(const std::size_t index) const;

                /// <summary>
                /// Return the specified original, advected seed point
                /// </summary>
                /// <param name="index">Index of the seed</param>
                /// <returns>Seed point</returns>
                tpf::geometry::point<floatp_t>& get_original_seed(const std::size_t index);

                /// <summary>
                /// Sort the particle traces for length and update the number of valid particles
                /// </summary>
                /// <param name="num_valid_particles">Previous number of valid particles</param>
                /// <returns>New number of valid particles</returns>
                virtual std::size_t sort_and_count(std::size_t num_valid_particles);

            protected:
                /// Advected seed positions
                std::vector<tpf::geometry::point<floatp_t>> original_seed;
            };
        }
    }
}

#include "tpf_streak_trace.inl"
