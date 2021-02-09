#include "tpf_path_trace.h"

#include "tpf_particle_seed.h"
#include "tpf_stream_trace.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/math/tpf_quaternion.h"

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
            inline path_trace<floatp_t>::path_trace(path_trace&& move) : stream_trace<floatp_t>(std::move(move))
            {
                this->original_particles = std::move(move.original_particles);
                this->rotations = std::move(move.rotations);
            }

            template <typename floatp_t>
            inline path_trace<floatp_t>::path_trace(particle_seed<floatp_t>&& move) : stream_trace<floatp_t>(std::move(move))
            {
                this->original_particles.resize(particle_seed<floatp_t>::seed.size());
                this->rotations.resize(particle_seed<floatp_t>::seed.size());

                for (std::size_t i = 0; i < particle_seed<floatp_t>::seed.size(); ++i)
                {
                    this->original_particles[i].push_back(particle_seed<floatp_t>::seed[i]);
                    this->rotations[i].resize(1);
                }
            }

            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& path_trace<floatp_t>::get_last_original_particle(const std::size_t index) const
            {
                return this->original_particles[index].back();
            }

            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& path_trace<floatp_t>::get_last_original_particle(const std::size_t index)
            {
                return this->original_particles[index].back();
            }

            template <typename floatp_t>
            inline const tpf::math::quaternion<floatp_t>& path_trace<floatp_t>::get_last_rotation(const std::size_t index) const
            {
                return this->rotations[index].back();
            }

            template <typename floatp_t>
            inline tpf::math::quaternion<floatp_t>& path_trace<floatp_t>::get_last_rotation(const std::size_t index)
            {
                return this->rotations[index].back();
            }

            template <typename floatp_t>
            inline void path_trace<floatp_t>::add_particle(const std::size_t index, const tpf::geometry::point<floatp_t>& particle,
                const tpf::geometry::point<floatp_t>& original_particle, const tpf::math::quaternion<floatp_t>& rotations)
            {
                stream_trace<floatp_t>::particles[index].push_back(particle);
                this->original_particles[index].push_back(original_particle);
                this->rotations[index].push_back(rotations);
            }

            template <typename floatp_t>
            inline std::size_t path_trace<floatp_t>::sort_and_count(const std::size_t num_valid_particles)
            {
                std::vector<geometry::point<floatp_t>> seed;
                std::vector<std::vector<geometry::point<floatp_t>>> particles;
                std::vector<std::vector<geometry::point<floatp_t>>> original_particles;
                std::vector<std::vector<math::quaternion<floatp_t>>> rotations;

                seed.reserve(particle_seed<floatp_t>::seed.size());
                particles.reserve(stream_trace<floatp_t>::particles.size());
                original_particles.reserve(this->original_particles.size());
                rotations.reserve(this->rotations.size());

                // Copy traces that are still valid
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (stream_trace<floatp_t>::validity[i])
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(stream_trace<floatp_t>::particles[i]);
                        original_particles.push_back(this->original_particles[i]);
                        rotations.push_back(this->rotations[i]);
                    }
                }

                const std::size_t new_num_valid_particles = particles.size();

                // Copy outdated traces
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (!stream_trace<floatp_t>::validity[i])
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(stream_trace<floatp_t>::particles[i]);
                        original_particles.push_back(this->original_particles[i]);
                        rotations.push_back(this->rotations[i]);
                    }
                }

                for (std::size_t i = num_valid_particles; i < stream_trace<floatp_t>::particles.size(); ++i)
                {
                    seed.push_back(particle_seed<floatp_t>::seed[i]);
                    particles.push_back(stream_trace<floatp_t>::particles[i]);
                    original_particles.push_back(this->original_particles[i]);
                    rotations.push_back(this->rotations[i]);
                }

                // Replace traces and seeds with reordered version
                std::swap(particle_seed<floatp_t>::seed, seed);
                std::swap(stream_trace<floatp_t>::particles, particles);
                std::swap(this->original_particles, original_particles);
                std::swap(this->rotations, rotations);

                return new_num_valid_particles;
            }
        }
    }
}
