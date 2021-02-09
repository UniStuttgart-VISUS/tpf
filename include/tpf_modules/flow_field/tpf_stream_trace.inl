#include "tpf_stream_trace.h"

#include "tpf_particle_seed.h"

#include "tpf/geometry/tpf_point.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            template <typename floatp_t>
            inline stream_trace<floatp_t>::stream_trace(stream_trace&& move) : particle_seed<floatp_t>(std::move(move))
            {
                this->particles = std::move(move.particles);

                this->validity.resize(this->particles.size());
                std::fill(this->validity.begin(), this->validity.end(), true);
            }

            template <typename floatp_t>
            inline stream_trace<floatp_t>::stream_trace(particle_seed<floatp_t>&& move) : particle_seed<floatp_t>(std::move(move))
            {
                this->particles.resize(particle_seed<floatp_t>::seed.size());
                this->validity.resize(this->particles.size());

                for (std::size_t i = 0; i < particle_seed<floatp_t>::seed.size(); ++i)
                {
                    this->particles[i].push_back(particle_seed<floatp_t>::seed[i]);
                    this->validity[i] = true;
                }
            }

            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& stream_trace<floatp_t>::get_last_particle(const std::size_t index) const
            {
                return this->particles[index].back();
            }
            
            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& stream_trace<floatp_t>::get_last_particle(const std::size_t index)
            {
                return this->particles[index].back();
            }

            template <typename floatp_t>
            inline void stream_trace<floatp_t>::add_particle(const std::size_t index, const tpf::geometry::point<floatp_t>& particle)
            {
                this->particles[index].push_back(particle);
            }

            template <typename floatp_t>
            inline std::size_t stream_trace<floatp_t>::get_num_lines() const
            {
                std::size_t num_lines = 0;

                for (const auto& trace : this->particles)
                {
                    num_lines += trace.size() - 1;
                }

                return num_lines;
            }

            template <typename floatp_t>
            inline void stream_trace<floatp_t>::invalidate_particle(const std::size_t index)
            {
                this->validity[index] = false;
            }

            template <typename floatp_t>
            inline const std::vector<std::vector<tpf::geometry::point<floatp_t>>>& stream_trace<floatp_t>::get_trace() const
            {
                return this->particles;
            }

            template <typename floatp_t>
            inline std::size_t stream_trace<floatp_t>::sort_and_count(const std::size_t num_valid_particles)
            {
                std::vector<geometry::point<floatp_t>> seed;
                std::vector<std::vector<geometry::point<floatp_t>>> particles;

                seed.reserve(particle_seed<floatp_t>::seed.size());
                particles.reserve(this->particles.size());

                // Copy traces that are still valid
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (this->validity[i])
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(this->particles[i]);
                    }
                }

                const std::size_t new_num_valid_particles = particles.size();

                // Copy outdated traces
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (!this->validity[i])
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(this->particles[i]);
                    }
                }

                for (std::size_t i = num_valid_particles; i < this->particles.size(); ++i)
                {
                    seed.push_back(particle_seed<floatp_t>::seed[i]);
                    particles.push_back(this->particles[i]);
                }

                // Replace traces and seeds with reordered version
                std::swap(this->particles, particles);
                std::swap(particle_seed<floatp_t>::seed, seed);

                return new_num_valid_particles;
            }
        }
    }
}
