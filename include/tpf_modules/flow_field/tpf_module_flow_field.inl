#include "tpf_module_flow_field.h"

#include "tpf_particle_seed.h"
#include "tpf_path_trace.h"
#include "tpf_streak_trace.h"
#include "tpf_stream_trace.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_line.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "Eigen/Dense"

#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t, typename point_t>
        inline std::size_t flow_field<float_t, point_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t, typename point_t>
        inline flow_field<float_t, point_t>::flow_field()
        { }

        template <typename float_t, typename point_t>
        inline std::string flow_field<float_t, point_t>::get_name() const
        {
            return std::string("Flow field");
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::set_algorithm_callbacks(flow_field_aux::request_frame_call_back<float_t, point_t>* next_time_frame_callback)
        {
            this->next_time_frame_callback = next_time_frame_callback;
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::set_algorithm_input(
            const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>& velocities,
            const policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, point_t>& global_velocity_parts,
            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_translation,
            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_rotation_axis,
            std::function<bool(const point_t&)> is_particle_valid, const data::polydata<float_t>& seed)
        {
            this->velocities = &velocities;
            this->global_velocity_parts = &global_velocity_parts;
            this->is_particle_valid = is_particle_valid;
            this->get_translation = get_translation;
            this->get_rotation_axis = get_rotation_axis;
            this->seed = &seed;
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::set_algorithm_output(tpf::data::polydata<float_t>& lines)
        {
            this->lines = &lines;
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::set_algorithm_parameters(const flow_field_aux::method_t method,
            const std::size_t num_advections, const float_t timestep)
        {
            this->method = method;
            this->num_advections = num_advections;
            this->timestep = timestep;
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::run_algorithm()
        {
            // Get seed
            const auto& seed_geometry = this->seed->get_geometry();

            std::vector<tpf::geometry::point<float_t>> seed_points;
            seed_points.reserve(seed_geometry.size());

            std::transform(seed_geometry.cbegin(), seed_geometry.cend(), std::back_inserter(seed_points), [](std::shared_ptr<geometry::geometric_object<float_t>> point)
                {
                    return static_cast<const tpf::geometry::point<float_t>&>(*point);
                });

            flow_field_aux::particle_seed<float_t> particles;
            particles.insert(seed_points.begin(), seed_points.end());

            // Compute integration lines
            if (this->method == flow_field_aux::method_t::STREAM)
            {
                compute_streamlines(std::move(particles));
            }
            else if (this->method == flow_field_aux::method_t::PATH)
            {
                compute_pathlines(std::move(particles));
            }
            else
            {
                compute_streaklines(std::move(particles));
            }
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::compute_streamlines(flow_field_aux::particle_seed<float_t>&& seed)
        {
            tpf::log::info_message(__tpf_info_message("Started streamline computation..."));

            // Create particle trace
            flow_field_aux::stream_trace<float_t> particles(std::move(seed));

            // Perform advection
            std::size_t advection = 0;
            std::size_t num_valid_particles = particles.get_size();

            for (; advection < this->num_advections && num_valid_particles > 0; ++advection)
            {
                // Duplicate and advect last particle
                /**
                    OpenMP information

                    write global:    particles
                **/
                #pragma omp parallel for schedule(static) default(none) shared(num_valid_particles, particles)
                for (long long index = 0; index < static_cast<long long>(num_valid_particles); ++index)
                {
                    const auto particle = particles.get_last_particle(index).get_vertex();

                    if (this->is_particle_valid(particle))
                    {
                        const auto velocity = this->velocities->interpolate(point_t(particle));
                        const auto global_velocity_part = this->global_velocity_parts->interpolate(point_t(particle));

                        const auto local_velocity = velocity - global_velocity_part;

                        if (!local_velocity.isZero())
                        {
                            particles.add_particle(index, tpf::geometry::point<float_t>(particle + local_velocity * this->timestep));
                        }
                    }
                }

                // Sort vector for number of particles and count still valid ones
                num_valid_particles = particles.sort_and_count(num_valid_particles, advection + 1);

                // Status information
                tpf::log::info_message(__tpf_info_message("Completed advection ", (advection + 1), " of ", this->num_advections));

#ifdef __tpf_debug
                tpf::log::info_message(__tpf_info_message("Number of valid particles: ", num_valid_particles));
#endif

                if (advection < this->num_advections - 1 && num_valid_particles == 0)
                {
                    tpf::log::info_message(__tpf_info_message("No valid particles for streamline integration. Completed prematurely after ", (advection + 1), " advections"));
                }
                else if (advection == this->num_advections - 1)
                {
                    tpf::log::info_message(__tpf_info_message("Completed streamline integration."));
                }
            }

            // Create line segments
            create_lines(particles, advection);
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::compute_pathlines(flow_field_aux::particle_seed<float_t>&& seed)
        {
            tpf::log::info_message(__tpf_info_message("Started pathline computation..."));

            // Get data and callback functions
            auto velocities = this->velocities;
            auto global_velocity_parts = this->global_velocity_parts;

            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_rotation_axis = this->get_rotation_axis;
            std::function<bool(const point_t&)> is_particle_valid = this->is_particle_valid;

            auto& next_time_frame_callback = *this->next_time_frame_callback;

            // Create particle trace
            flow_field_aux::path_trace<float_t> particles(std::move(seed));

            // Perform advection
            std::size_t advection = 0;
            std::size_t num_valid_particles = particles.get_size();
            bool has_available_timestep = true;

            for (; advection < this->num_advections && num_valid_particles > 0 && has_available_timestep; ++advection)
            {
                // Duplicate and advect last particle
                /**
                    OpenMP information

                    write global:    particles
                **/
                #pragma omp parallel for schedule(static) default(none) shared(num_valid_particles, particles, velocities, \
                    is_particle_valid, get_rotation_axis, global_velocity_parts)
                for (long long index = 0; index < static_cast<long long>(num_valid_particles); ++index)
                {
                    const auto particle = particles.get_last_particle(index).get_vertex();
                    const auto original_particle = particles.get_last_original_particle(index).get_vertex();

                    if (is_particle_valid(original_particle))
                    {
                        const auto velocity = velocities->interpolate(point_t(original_particle));
                        const auto global_velocity_part = global_velocity_parts->interpolate(point_t(original_particle));

                        Eigen::Matrix<float_t, 3, 1> local_velocity = velocity - global_velocity_part;

                        if (!local_velocity.isZero())
                        {
                            // Calculate rotation
                            const auto previous_rotation = particles.get_last_rotation(index);
                            const auto axis = get_rotation_axis(original_particle);

                            tpf::math::quaternion<float_t> quaternion;
                            quaternion.from_axis(axis);
                            quaternion.invert();

                            const auto new_rotation = quaternion * previous_rotation;

                            // Calculate droplet-local velocity
                            local_velocity = previous_rotation * local_velocity;

                            particles.add_particle(index, tpf::geometry::point<float_t>(particle + local_velocity * this->timestep),
                                tpf::geometry::point<float_t>(original_particle + velocity * this->timestep), new_rotation);
                        }
                    }
                }

                // Sort vector for number of particles and count still valid ones
                num_valid_particles = particles.sort_and_count(num_valid_particles, advection + 1);

                // Status information
                tpf::log::info_message(__tpf_info_message("Completed advection ", (advection + 1), " of ", this->num_advections));

#ifdef __tpf_debug
                tpf::log::info_message(__tpf_info_message("Number of valid particles: ", num_valid_particles));
#endif

                if (advection < this->num_advections - 1 && num_valid_particles == 0)
                {
                    tpf::log::info_message(__tpf_info_message("No valid particles for pathline integration. Completed prematurely after ", (advection + 1), " advections"));
                }
                else if (advection == this->num_advections - 1)
                {
                    tpf::log::info_message(__tpf_info_message("Completed pathline integration."));
                }

                // Try to load the next time step
                try
                {
                    if (advection < this->num_advections - 1 && num_valid_particles > 0)
                    {
                        // Load next time step's data
                        try
                        {
                            // Get data for next time step: [Time step delta, velocities]
                            auto next_data = next_time_frame_callback();

                            this->timestep = std::get<0>(next_data) == static_cast<float_t>(0.0) ? this->timestep : std::get<0>(next_data);
                            velocities = std::get<1>(next_data);
                            global_velocity_parts = std::get<2>(next_data);
                            get_rotation_axis = std::get<4>(next_data);
                            is_particle_valid = std::get<5>(next_data);
                        }
                        catch (const std::exception&)
                        {
                            has_available_timestep = false;
                        }
                    }

                    if (advection == this->num_advections - 1 || !has_available_timestep || num_valid_particles == 0)
                    {
                        // Reset data-providing modules
                        next_time_frame_callback.reset();
                    }
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error occured while loading next time frame."));
                }
            }

            // Create line segments
            create_lines(particles, advection);
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::compute_streaklines(flow_field_aux::particle_seed<float_t>&& seed)
        {
            tpf::log::info_message(__tpf_info_message("Started streakline computation..."));

            // Get data for whole extent
            auto velocities = this->velocities;
            auto global_velocity_parts = this->global_velocity_parts;

            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_translation = this->get_translation;
            std::function<Eigen::Matrix<float_t, 3, 1>(const point_t&)> get_rotation_axis = this->get_rotation_axis;
            std::function<bool(const point_t&)> is_particle_valid = this->is_particle_valid;

            auto& next_time_frame_callback = *this->next_time_frame_callback;

            // Create particle trace
            flow_field_aux::streak_trace<float_t> particles(std::move(seed));

            // Perform advection
            std::size_t advection = 0;
            std::size_t num_valid_particles = particles.get_size();
            bool has_available_timestep = true;

            for (; advection < this->num_advections && num_valid_particles > 0 && has_available_timestep; ++advection)
            {
                // Advect all particles and add new seed particle
                /**
                    OpenMP information

                    write global:    particles
                **/
                #pragma omp parallel for schedule(static) default(none) shared(num_valid_particles, particles, velocities \
                    is_particle_valid, get_rotation_axis, get_translation, global_velocity_parts)
                for (long long index = 0; index < static_cast<long long>(num_valid_particles); ++index)
                {
                    for (std::size_t particle_index = 0; particle_index < particles.get_num_particles(index); ++particle_index)
                    {
                        const auto particle = particles.get_particle(index, particle_index).get_vertex();
                        const auto original_particle = particles.get_original_particle(index, particle_index).get_vertex();

                        if (is_particle_valid(original_particle))
                        {
                            const auto velocity = velocities->interpolate(point_t(original_particle));
                            const auto global_velocity_part = global_velocity_parts->interpolate(point_t(original_particle));

                            Eigen::Matrix<float_t, 3, 1> local_velocity = velocity - global_velocity_part;

                            if (!local_velocity.isZero())
                            {
                                // Calculate rotation
                                const auto previous_rotation = particles.get_rotation(index, particle_index);
                                const auto axis = get_rotation_axis(original_particle);

                                tpf::math::quaternion<float_t> quaternion;
                                quaternion.from_axis(axis);
                                quaternion.invert();

                                const auto new_rotation = quaternion * previous_rotation;

                                // Calculate droplet-local velocity
                                local_velocity = previous_rotation * local_velocity;

                                // Advect particles
                                particles.get_particle(index, particle_index) = tpf::geometry::point<float_t>(particle + local_velocity * this->timestep);
                                particles.get_original_particle(index, particle_index) = tpf::geometry::point<float_t>(original_particle + velocity * this->timestep);
                                particles.get_rotation(index, particle_index) = new_rotation;
                            }
                        }
                    }

                    // Interpolate velocity and advect the original seed, while keeping the local seed unchanged
                    auto& original_seed = particles.get_original_seed(index);

                    if (is_particle_valid(original_seed))
                    {
                        // Advect the seed using the global velocity part
                        const auto translation = get_translation(original_seed);
                        const auto axis = get_rotation_axis(original_seed);

                        tpf::math::quaternion<float_t> quaternion;
                        quaternion.from_axis(axis);
                        quaternion.invert();

                        original_seed = tpf::geometry::point<float_t>(original_seed.get_vertex()
                            + (((quaternion * original_seed.get_vertex()) - original_seed.get_vertex()) + translation) * this->timestep);

                        // Insert new particle at the seed
                        particles.add_particle(index, particles.get_seed(index), original_seed, tpf::math::quaternion<float_t>());
                    }
                }

                // Sort vector for number of particles and count still valid ones
                num_valid_particles = particles.sort_and_count(num_valid_particles, advection + 1);

                // Status information
                tpf::log::info_message(__tpf_info_message("Completed advection ", (advection + 1), " of ", this->num_advections));

#ifdef __tpf_debug
                tpf::log::info_message(__tpf_info_message("Number of valid particles: ", num_valid_particles));
#endif

                if (advection < this->num_advections - 1 && num_valid_particles == 0)
                {
                    tpf::log::info_message(__tpf_info_message("No valid particles for streakline integration. Completed prematurely after ", (advection + 1), " advections"));
                }
                else if (advection == this->num_advections - 1)
                {
                    tpf::log::info_message(__tpf_info_message("Completed streakline integration."));
                }

                // Try to load the next time step
                try
                {
                    if (advection < this->num_advections - 1 && num_valid_particles > 0)
                    {
                        // Load next time step's data
                        try
                        {
                            // Get data for next time step: [Time step delta, velocities]
                            auto next_data = next_time_frame_callback();

                            this->timestep = std::get<0>(next_data) == static_cast<float_t>(0.0) ? this->timestep : std::get<0>(next_data);
                            velocities = std::get<1>(next_data);
                            global_velocity_parts = std::get<2>(next_data);
                            get_translation = std::get<3>(next_data);
                            get_rotation_axis = std::get<4>(next_data);
                            is_particle_valid = std::get<5>(next_data);
                        }
                        catch (const std::exception&)
                        {
                            has_available_timestep = false;
                        }
                    }

                    if (advection == this->num_advections - 1 || !has_available_timestep || num_valid_particles == 0)
                    {
                        // Reset data-providing modules
                        next_time_frame_callback.reset();
                    }
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error occured while loading next time frame."));
                }
            }

            // Create line segments
            create_lines(particles, advection, true);
        }

        template <typename float_t, typename point_t>
        inline void flow_field<float_t, point_t>::create_lines(const flow_field_aux::stream_trace<float_t>& particles, const std::size_t num_advections, const bool inverse)
        {
            // Create array
            auto id_advection_ref = std::make_shared<tpf::data::array<std::size_t, 1>>("ID (Advection)", particles.get_num_lines());
            auto id_distribution_ref = std::make_shared<tpf::data::array<std::size_t, 1>>("ID (Distribution)", particles.get_num_lines());

            auto& id_advection = *id_advection_ref;
            auto& id_distribution = *id_distribution_ref;

            std::size_t index = 0;

            // Create line segments
            for (const auto& trace : particles.get_trace())
            {
                const auto distribution_step = static_cast<float_t>(num_advections) / static_cast<float_t>(trace.size());

                if (!inverse)
                {
                    for (std::size_t i = 0; i < trace.size() - 1; ++i)
                    {
                        this->lines->insert(std::make_shared<tpf::geometry::line<float_t>>(trace[i], trace[i + 1]));

                        id_advection(index) = i;
                        id_distribution(index) = static_cast<std::size_t>(i * distribution_step);

                        ++index;
                    }
                }
                else
                {
                    for (std::size_t i = trace.size() - 1; i > 0; --i)
                    {
                        this->lines->insert(std::make_shared<tpf::geometry::line<float_t>>(trace[i], trace[i - 1]));

                        id_advection(index) = i;
                        id_distribution(index) = static_cast<std::size_t>(i * distribution_step);

                        ++index;
                    }
                }
            }

            // Store directional line information (IDs)
            this->lines->add(id_advection_ref, tpf::data::topology_t::CELL_DATA);
            this->lines->add(id_distribution_ref, tpf::data::topology_t::CELL_DATA);
        }
    }
}
