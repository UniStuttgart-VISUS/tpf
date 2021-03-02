#pragma once

#include "tpf_particle_seed.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/policies/tpf_interpolatable.h"

#include "tpf/stdext/tpf_bool.h"

#include <string>
#include <utility>
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
            class stream_trace : public particle_seed<floatp_t>
            {
            public:
                /// <summary>
                /// Disallow use of copy constructor
                /// </summary>
                stream_trace(const stream_trace&) = delete;

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source object</param>
                stream_trace(stream_trace&& move);

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source particle seed object to create this from</param>
                /// <param name="property_names">Name of attached property arrays</param>
                /// <param name="properties">Properties interpolated at the respective seed</param>
                stream_trace(particle_seed<floatp_t>&& move, std::vector<std::string>&& property_names,
                    std::vector<std::vector<std::vector<double>>>&& properties);

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
                /// Add particle to a given trace
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                /// <param name="particle">New particle</param>
                /// <param name="properties">Properties interpolated at the particle</param>
                void add_particle(std::size_t index, const tpf::geometry::point<floatp_t>& particle,
                    std::vector<std::vector<double>>&& properties);

                /// <summary>
                /// Return the number of lines, given by the number of traces and their number of particles
                /// </summary>
                /// <returns>Number of lines</returns>
                std::size_t get_num_lines() const;

                /// <summary>
                /// Return the number of (advected) particles
                /// </summary>
                /// <returns>Number of particles</returns>
                std::size_t get_num_particles() const;

                /// <summary>
                /// Invalidate a particle (trace)
                /// </summary>
                /// <param name="index">Index of the seed/trace</param>
                void invalidate_particle(std::size_t index);

                /// <summary>
                /// Sort the particle traces for length and update the number of valid particles
                /// </summary>
                /// <param name="num_valid_particles">Previous number of valid particles</param>
                /// <returns>New number of valid particles</returns>
                virtual std::size_t sort_and_count(std::size_t num_valid_particles);

                /// <summary>
                /// Return the traces
                /// </summary>
                /// <returns>Traces</returns>
                const std::vector<std::vector<tpf::geometry::point<floatp_t>>>& get_trace() const;

                /// <summary>
                /// Return particle properties for all traces
                /// </summary>
                /// <returns>Particle properties</returns>
                const std::vector<std::vector<std::vector<std::vector<double>>>>& get_particle_properties() const;
                const std::vector<std::string>& get_particle_property_names() const;

            protected:
                /// Particles
                std::vector<std::vector<tpf::geometry::point<floatp_t>>> particles;

                /// Validity of particles
                std::vector<tpf::bool_t> validity;

                /// Store properties at particles
                std::vector<std::vector<std::vector<std::vector<double>>>> particle_properties;
                std::vector<std::string> property_names;
            };
        }
    }
}

#include "tpf_stream_trace.inl"
