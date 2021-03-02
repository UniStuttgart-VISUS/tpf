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
                this->validity = std::move(move.validity);
                this->seed_properties = std::move(move.seed_properties);
                this->particle_properties = std::move(move.particle_properties);
                this->property_names = std::move(move.property_names);
            }

            template <typename floatp_t>
            inline stream_trace<floatp_t>::stream_trace(particle_seed<floatp_t>&& move, std::vector<std::string>&& property_names,
                std::vector<std::vector<std::vector<double>>>&& properties) : particle_seed<floatp_t>(std::move(move))
            {
                this->particles.resize(particle_seed<floatp_t>::seed.size());
                this->validity.resize(this->particles.size());

                this->particle_properties.resize(properties.size());
                this->property_names = std::move(property_names);

                std::for_each(this->particle_properties.begin(), this->particle_properties.end(),
                    [this](std::vector<std::vector<std::vector<double>>>& property) { property.resize(this->particles.size()); });

                for (std::size_t i = 0; i < particle_seed<floatp_t>::seed.size(); ++i)
                {
                    this->particles[i].push_back(particle_seed<floatp_t>::seed[i]);
                    this->validity[i] = true;

                    for (std::size_t p = 0; p < properties.size(); ++p)
                    {
                        this->particle_properties[p][i].push_back(std::move(properties[p][i]));
                    }
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
            inline void stream_trace<floatp_t>::add_particle(const std::size_t index, const tpf::geometry::point<floatp_t>& particle,
                std::vector<std::vector<double>>&& properties)
            {
                this->particles[index].push_back(particle);

                for (std::size_t p = 0; p < this->particle_properties.size(); ++p)
                {
                    this->particle_properties[p][index].push_back(std::move(properties[p]));
                }
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
            inline std::size_t stream_trace<floatp_t>::get_num_particles() const
            {
                std::size_t num_particles = 0;

                for (const auto& trace : this->particles)
                {
                    num_particles += trace.size();
                }

                return num_particles;
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
            inline const std::vector<std::vector<std::vector<std::vector<double>>>>& stream_trace<floatp_t>::get_particle_properties() const
            {
                return this->particle_properties;
            }

            template <typename floatp_t>
            inline const std::vector<std::string>& stream_trace<floatp_t>::get_particle_property_names() const
            {
                return this->property_names;
            }

            template <typename floatp_t>
            inline std::size_t stream_trace<floatp_t>::sort_and_count(const std::size_t num_valid_particles)
            {
                std::vector<geometry::point<floatp_t>> seed;
                std::vector<std::vector<geometry::point<floatp_t>>> particles;
                std::vector<std::vector<std::vector<std::vector<double>>>> particle_properties;

                seed.reserve(particle_seed<floatp_t>::seed.size());
                particles.reserve(this->particles.size());

                particle_properties.resize(this->particle_properties.size());

                std::for_each(particle_properties.begin(), particle_properties.end(),
                    [this](std::vector<std::vector<std::vector<double>>>& property) { property.reserve(this->particles.size()); });

                // Copy traces that are still valid
                for (std::size_t i = 0; i < num_valid_particles; ++i)
                {
                    if (this->validity[i])
                    {
                        seed.push_back(particle_seed<floatp_t>::seed[i]);
                        particles.push_back(this->particles[i]);

                        for (std::size_t p = 0; p < particle_properties.size(); ++p)
                        {
                            particle_properties[p].push_back(this->particle_properties[p][i]);
                        }
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

                        for (std::size_t p = 0; p < particle_properties.size(); ++p)
                        {
                            particle_properties[p].push_back(this->particle_properties[p][i]);
                        }
                    }
                }

                for (std::size_t i = num_valid_particles; i < this->particles.size(); ++i)
                {
                    seed.push_back(particle_seed<floatp_t>::seed[i]);
                    particles.push_back(this->particles[i]);

                    for (std::size_t p = 0; p < particle_properties.size(); ++p)
                    {
                        particle_properties[p].push_back(this->particle_properties[p][i]);
                    }
                }

                // Replace traces and seeds with reordered version
                std::swap(particle_seed<floatp_t>::seed, seed);
                std::swap(this->particles, particles);
                std::swap(this->particle_properties, particle_properties);

                return new_num_valid_particles;
            }
        }
    }
}
