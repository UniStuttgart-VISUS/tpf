#pragma once

#include "tpf_particle_seed.h"

#include "tpf/geometry/tpf_point.h"

#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            /// <summary>
            /// Class for storing particle traces and seeds
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class simple_trace : public particle_seed<floatp_t>
            {
            public:
                /// <summary>
                /// Disallow use of copy constructor
                /// </summary>
                simple_trace(const simple_trace&) = delete;

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source object</param>
                simple_trace(simple_trace&& move);

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source particle seed object to create this from</param>
                simple_trace(particle_seed<floatp_t>&& move);

                /// <summary>
                /// Return last particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last particle point</returns>
                const tpf::geometry::point<floatp_t>& get_last_particle(std::size_t index) const;

                /// <summary>
                /// Return last particle for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Last particle point</returns>
                tpf::geometry::point<floatp_t>& get_last_particle(std::size_t index);

                /// <summary>
                /// Return particles for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Particle points</returns>
                const std::vector<tpf::geometry::point<floatp_t>>& get_particles(std::size_t index) const;

                /// <summary>
                /// Return particles for a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <returns>Particle points</returns>
                std::vector<tpf::geometry::point<floatp_t>>& get_particles(std::size_t index);

                /// <summary>
                /// Add particle to a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle">New particle</param>
                void add_particle(std::size_t index, const tpf::geometry::point<floatp_t>& particle);

                /// <summary>
                /// Return the number of lines, given by the number of traces and their number of particles
                /// </summary>
                /// <returns>Number of lines</returns>
                std::size_t get_num_lines() const;

                /// <summary>
                /// Sort the particle traces for length and update the number of valid particles
                /// </summary>
                /// <param name="num_valid_particles">Previous number of valid particles</param>
                /// <param name="num_advections">Number of advection step</param>
                /// <returns>New number of valid particles</returns>
                virtual std::size_t sort_and_count(std::size_t num_valid_particles, std::size_t num_advections);

                /// <summary>
                /// Return the traces
                /// </summary>
                /// <returns>Traces</returns>
                const std::vector<std::vector<tpf::geometry::point<floatp_t>>>& get_trace() const;

            protected:
                /// Particles
                std::vector<std::vector<tpf::geometry::point<floatp_t>>> particles;
            };
        }
    }
}

#include "tpf_trace.inl"
