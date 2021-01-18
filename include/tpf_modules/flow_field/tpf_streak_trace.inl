#include "tpf_streak_trace.h"

#include "tpf_particle_seed.h"
#include "tpf_path_trace.h"

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
            inline streak_trace<floatp_t>::streak_trace(streak_trace&& move) : path_trace<floatp_t>(std::move(move))
            {
                this->original_seed = std::move(move.original_seed);
            }

            template <typename floatp_t>
            inline streak_trace<floatp_t>::streak_trace(particle_seed<floatp_t>&& move) : path_trace<floatp_t>(std::move(move))
            {
                this->original_seed = particle_seed<floatp_t>::seed;
            }

            template <typename floatp_t>
            inline std::size_t streak_trace<floatp_t>::get_num_particles(const std::size_t index) const
            {
                return stream_trace<floatp_t>::particles[index].size();
            }

            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_particle(const std::size_t index, const std::size_t particle_index) const
            {
                return stream_trace<floatp_t>::particles[index][particle_index];
            }

            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_particle(const std::size_t index, const std::size_t particle_index)
            {
                return stream_trace<floatp_t>::particles[index][particle_index];
            }

            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_original_particle(const std::size_t index, const std::size_t particle_index) const
            {
                return path_trace<floatp_t>::original_particles[index][particle_index];
            }

            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_original_particle(const std::size_t index, const std::size_t particle_index)
            {
                return path_trace<floatp_t>::original_particles[index][particle_index];
            }

            template <typename floatp_t>
            inline const tpf::math::quaternion<floatp_t>& streak_trace<floatp_t>::get_rotation(const std::size_t index, const std::size_t particle_index) const
            {
                return path_trace<floatp_t>::rotations[index][particle_index];
            }

            template <typename floatp_t>
            inline tpf::math::quaternion<floatp_t>& streak_trace<floatp_t>::get_rotation(const std::size_t index, const std::size_t particle_index)
            {
                return path_trace<floatp_t>::rotations[index][particle_index];
            }

            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_original_seed(const std::size_t index) const
            {
                return this->original_seed[index];
            }

            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& streak_trace<floatp_t>::get_original_seed(const std::size_t index)
            {
                return this->original_seed[index];
            }

            template <typename floatp_t>
            inline std::size_t streak_trace<floatp_t>::sort_and_count(const std::size_t num_valid_particles, const std::size_t num_advections)
            {
                std::vector<geometry::point<floatp_t>> seed;
                std::vector<std::vector<geometry::point<floatp_t>>> particles;
                std::vector<std::vector<geometry::point<floatp_t>>> original_particles;
                std::vector<std::vector<math::quaternion<floatp_t>>> rotations;
                std::vector<geometry::point<floatp_t>> original_seed;

                seed.reserve(particle_seed<floatp_t>::seed.size());
                particles.reserve(stream_trace<floatp_t>::particles.size());
                original_particles.reserve(path_trace<floatp_t>::original_particles.size());
                rotations.reserve(path_trace<floatp_t>::rotations.size());
                original_seed.reserve(this->original_seed.size());

                // Copy traces that are still valid
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (stream_trace<floatp_t>::particles[i].size() == num_advections + 1)
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(stream_trace<floatp_t>::particles[i]);
                        original_particles.push_back(path_trace<floatp_t>::original_particles[i]);
                        rotations.push_back(path_trace<floatp_t>::rotations[i]);
                        original_seed.push_back(this->original_seed[i]);
                    }
                }

                const std::size_t new_num_valid_particles = particles.size();

                // Copy outdated traces
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (stream_trace<floatp_t>::particles[i].size() != num_advections + 1)
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(stream_trace<floatp_t>::particles[i]);
                        original_particles.push_back(path_trace<floatp_t>::original_particles[i]);
                        rotations.push_back(path_trace<floatp_t>::rotations[i]);
                        original_seed.push_back(this->original_seed[i]);
                    }
                }

                for (std::size_t i = num_valid_particles; i < stream_trace<floatp_t>::particles.size(); ++i)
                {
                    seed.push_back(particle_seed<floatp_t>::seed[i]);
                    particles.push_back(stream_trace<floatp_t>::particles[i]);
                    original_particles.push_back(path_trace<floatp_t>::original_particles[i]);
                    rotations.push_back(path_trace<floatp_t>::rotations[i]);
                    original_seed.push_back(this->original_seed[i]);
                }

                // Replace traces and seeds with reordered version
                std::swap(particle_seed<floatp_t>::seed, seed);
                std::swap(stream_trace<floatp_t>::particles, particles);
                std::swap(path_trace<floatp_t>::original_particles, original_particles);
                std::swap(path_trace<floatp_t>::rotations, rotations);
                std::swap(this->original_seed, original_seed);

                return new_num_valid_particles;
            }
        }
    }
}
