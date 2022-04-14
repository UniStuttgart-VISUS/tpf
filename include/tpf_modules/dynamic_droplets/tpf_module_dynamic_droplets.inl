#include "tpf_module_dynamic_droplets.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_line.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_quaternion.h"
#include "tpf/math/tpf_vector.h"

#include "tpf/stdext/tpf_comparator.h"

#include "Eigen/Dense"

#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t dynamic_droplets<float_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t>
        inline dynamic_droplets<float_t>::dynamic_droplets()
        { }

        template <typename float_t>
        inline std::string dynamic_droplets<float_t>::get_name() const
        {
            return std::string("Dynamic droplets");
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_callbacks(dynamic_droplets_aux::request_frame_call_back<float_t>* next_time_frame_callback)
        {
            this->next_time_frame_callback = next_time_frame_callback;
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_input(const data::polydata<float_t>& droplets,
            const data::grid<long long, float_t, 3, 1>& droplet_ids)
        {
            this->droplets = &droplets;
            this->droplet_ids = &droplet_ids;
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_output(data::polydata<float_t>& tracks,
            data::polydata<float_t>& summary, data::polydata<float_t>& paths, data::polydata<float_t>& axes,
            data::polydata<float_t>& ribbons, data::polydata<float_t>& rotation_paths, data::polydata<float_t>& coordinate_axes)
        {
            this->tracks = &tracks;
            this->summary = &summary;
            this->paths = &paths;
            this->axes = &axes;
            this->ribbons = &ribbons;
            this->rotation_paths = &rotation_paths;
            this->coordinate_axes = &coordinate_axes;
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_parameters(const std::size_t num_timesteps,
            const float_t timestep, const bool static_frame_of_reference, const float_t ribbon_scale,
            const float_t ribbon_thickness, const dynamic_droplets_aux::ribbon_duplication_t ribbon_duplication,
            const bool translation_duplication, const bool fix_axis_size, const float_t axis_scale, const bool axis_translation,
            std::string translation_name, std::string rotation_name, std::string radius_name)
        {
            this->num_timesteps = num_timesteps;
            this->timestep = timestep;

            this->static_frame_of_reference = static_frame_of_reference;
            this->ribbon_scale = ribbon_scale;
            this->ribbon_thickness = ribbon_thickness;
            this->ribbon_duplication = ribbon_duplication;
            this->translation_duplication = translation_duplication;
            this->fix_axis_size = fix_axis_size;
            this->axis_scale = axis_scale;
            this->axis_translation = axis_translation;

            this->translation_name = translation_name;
            this->rotation_name = rotation_name;
            this->radius_name = radius_name;
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::run_algorithm()
        {
            // Check if there is anything to do
            if (this->num_timesteps <= 1 || this->droplets->get_num_objects() == 0)
            {
                log::info_message(__tpf_info_message("No data. Nothing to do."));

                return;
            }

            // Track droplets over time
            const auto droplets_over_time = track_droplets();

            // Output tracking information
            this->tracks->add(std::make_shared<data::array<std::size_t>>("Topology information"), data::topology_t::CELL_DATA);

            for (const auto& droplet_track : droplets_over_time.first)
            {
                const auto topology = droplet_track.second;

                for (std::size_t i = 0; i < droplet_track.first.size() - 1; ++i)
                {
                    this->tracks->insert(std::make_shared<geometry::line<float_t>>(
                        droplet_track.first[i].position, droplet_track.first[i + 1].position),
                        data::polydata_element<std::size_t, 1>("Topology information", data::topology_t::CELL_DATA,
                            { static_cast<std::size_t>(topology) }));
                }
            }

            // Create output
            create_summary(droplets_over_time.first, droplets_over_time.second);
            create_paths(droplets_over_time.first, droplets_over_time.second);
            create_axes(droplets_over_time.first, droplets_over_time.second);
            create_ribbons(droplets_over_time.first, droplets_over_time.second);
            create_rotation_paths(droplets_over_time.first, droplets_over_time.second);
            create_coordinate_axes(droplets_over_time.first, droplets_over_time.second);
        }

        template <typename float_t>
        inline std::pair<std::vector<typename dynamic_droplets<float_t>::droplet_trace_t>, std::vector<float_t>>
            dynamic_droplets<float_t>::track_droplets() const
        {
            // Collect droplets for all time steps
            auto timestep_delta = this->timestep;
            auto droplets = *this->droplets;
            auto droplet_ids = *this->droplet_ids;

            std::vector<float_t> timesteps;
            timesteps.reserve(this->num_timesteps);
            timesteps.push_back(timestep_delta);

            // Map droplets from different time steps, and record topological changes
            std::vector<droplet_trace_t> droplets_over_time(this->droplets->get_num_objects());

            std::map<std::size_t, std::size_t> map_to_original;

            auto position = this->droplets->get_geometry();
            auto translation = this->droplets->template get_point_data_as<float_t, 3>(this->translation_name);
            auto rotation = this->droplets->template get_point_data_as<float_t, 3>(this->rotation_name);
            auto radius = this->droplets->template get_point_data_as<float_t, 1>(this->radius_name);

            for (std::size_t i = 0; i < droplets_over_time.size(); ++i)
            {
                droplets_over_time[i].first.push_back(droplet_t
                    { std::dynamic_pointer_cast<geometry::point<float_t>>(position[i])->get_vertex(),
                    (*translation)(i), (*rotation)(i), (*radius)(i) });

                droplets_over_time[i].second = topology_t::none;

                map_to_original[i] = i;
            }

            for (std::size_t timestep = 1; timestep <= this->num_timesteps; ++timestep)
            {
                // Try to load the next time step
                bool has_available_timestep = true;

                try
                {
                    if (timestep < this->num_timesteps)
                    {
                        // Load next time step's data
                        try
                        {
                            auto next_data = (*this->next_time_frame_callback)();

                            // Get complete data: [Time step delta, droplets]
                            auto& next_timestep_delta = std::get<0>(next_data);
                            auto& next_droplets = std::get<1>(next_data);
                            auto& next_droplet_ids = std::get<2>(next_data);

                            auto next_position = next_droplets.get_geometry();
                            auto next_translation = next_droplets.template get_point_data_as<float_t, 3>(this->translation_name);
                            auto next_rotation = next_droplets.template get_point_data_as<float_t, 3>(this->rotation_name);
                            auto next_radius = next_droplets.template get_point_data_as<float_t, 1>(this->radius_name);

                            timesteps.push_back(next_timestep_delta);

                            // Find corresponding future droplets by forward advection -> detects collision
                            std::map<std::size_t, std::vector<std::size_t>> forw_adv_backw_mapping;

                            {
                                // Create map from future to current droplet
                                for (std::size_t i = 0; i < droplets.get_num_objects(); ++i)
                                {
                                    if (map_to_original.find(i) != map_to_original.end() &&
                                        droplets_over_time[map_to_original.at(i)].second == topology_t::none)
                                    {
                                        const auto advected_position = std::dynamic_pointer_cast<geometry::point<float_t>>(position[i])->get_vertex()
                                            + timestep_delta * (*translation)(i);

                                        const auto cell = next_droplet_ids.find_cell(advected_position);

                                        if (cell && next_droplet_ids(*cell) >= 0)
                                        {
                                            const auto droplet_id = next_droplet_ids(*cell);

                                            forw_adv_backw_mapping[droplet_id].push_back(i);
                                        }
                                    }
                                }
                            }

                            // Find corresponding current droplets by backward advection -> detects breakup
                            std::map<std::size_t, std::vector<std::size_t>> backw_adv_forw_mapping;

                            {
                                // Create map from current to future droplet
                                for (std::size_t i = 0; i < next_droplets.get_num_objects(); ++i)
                                {
                                    const auto advected_position = std::dynamic_pointer_cast<geometry::point<float_t>>(next_position[i])->get_vertex()
                                        - next_timestep_delta * (*next_translation)(i);

                                    const auto cell = droplet_ids.find_cell(advected_position);

                                    if (cell && droplet_ids(*cell) >= 0)
                                    {
                                        const auto droplet_id = droplet_ids(*cell);

                                        if (map_to_original.find(droplet_id) != map_to_original.end())
                                        {
                                            backw_adv_forw_mapping[droplet_id].push_back(i);
                                        }
                                    }
                                }
                            }

                            // Classify topology changes: 0 - none, 1 - breakup, 2 - collision
                            std::map<std::size_t, std::size_t> new_map_to_original;

                            for (const auto& breakup_map : backw_adv_forw_mapping)
                            {
                                // Breakup
                                if (breakup_map.second.size() > 1)
                                {
                                    droplets_over_time[map_to_original.at(breakup_map.first)].second = topology_t::breakup;
                                }
                            }

                            for (const auto& collision_map : forw_adv_backw_mapping)
                            {
                                // Collision
                                if (collision_map.second.size() > 1)
                                {
                                    for (std::size_t i = 0; i < collision_map.second.size(); ++i)
                                    {
                                        droplets_over_time[map_to_original.at(collision_map.second[i])].second = topology_t::collision;
                                    }
                                }
                                else
                                {
                                    // No topological change: add droplet
                                    const auto original_droplet = map_to_original.at(collision_map.second.front());
                                    new_map_to_original[collision_map.first] = original_droplet;

                                    if (!this->static_frame_of_reference)
                                    {
                                        droplets_over_time[original_droplet].first.push_back(droplet_t{
                                            std::dynamic_pointer_cast<geometry::point<float_t>>(next_position[collision_map.first])->get_vertex(),
                                            (*next_translation)(collision_map.first),
                                            (*next_rotation)(collision_map.first),
                                            (*next_radius)(collision_map.first)
                                            });
                                    }
                                    else
                                    {
                                        // Add pseudo droplet as if the frame of reference was static
                                        droplets_over_time[original_droplet].first.push_back(droplet_t{
                                            droplets_over_time[original_droplet].first.front().position
                                                + timestep * timesteps.front() * droplets_over_time[original_droplet].first.front().translation,
                                            droplets_over_time[original_droplet].first.front().translation,
                                            droplets_over_time[original_droplet].first.front().rotation,
                                            droplets_over_time[original_droplet].first.front().radius
                                            });
                                    }
                                }
                            }

                            // Prepare for next time step
                            position = next_position;
                            translation = next_translation;
                            rotation = next_rotation;
                            radius = next_radius;

                            std::swap(map_to_original, new_map_to_original);

                            std::swap(timestep_delta, next_timestep_delta);
                            std::swap(droplets, next_droplets);
                            std::swap(droplet_ids, next_droplet_ids);
                        }
                        catch (const std::exception&)
                        {
                            has_available_timestep = false;
                        }
                    }

                    if (timestep == this->num_timesteps || !has_available_timestep)
                    {
                        // Reset data-providing modules
                        this->next_time_frame_callback->reset();
                    }
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error occured while loading next time frame."));
                }

                log::info_message(__tpf_info_message("Completed time step ", timestep, " of ", this->num_timesteps, "."));

#ifdef __tpf_debug
                log::info_message(__tpf_info_message("  Number of droplets: ", droplets.get_geometry().size()));
                log::info_message(__tpf_info_message("  Time step size: ", timestep_delta));
#endif

                if (timestep == this->num_timesteps)
                {
                    log::info_message(__tpf_info_message("Completed loading all time steps."));
                }
            }

            return std::make_pair(droplets_over_time, timesteps);
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_summary(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Create position with displacement information
            const std::string data_name("Displacement");

            this->summary->add(std::make_shared<data::array<float_t, 3>>(data_name), data::topology_t::POINT_DATA);

            for (const auto& droplet : all_droplets)
            {
                const auto position = droplet.first.front().position;
                const auto radius = droplet.first.front().radius;

                const Eigen::Matrix<float_t, 3, 1> direction = droplet.first.front().translation.normalized();
                const Eigen::Matrix<float_t, 3, 1> start_position = position + radius * direction;

                auto droplet_position = start_position;

                for (std::size_t i = 0; i < droplet.first.size(); ++i)
                {
                    droplet_position += timesteps[i] * droplet.first[i].translation;
                }

                this->summary->insert(std::make_shared<geometry::point<float_t>>(start_position),
                    data::polydata_element<float_t, 3>(data_name, data::topology_t::POINT_DATA,
                        { Eigen::Matrix<float_t, 3, 1>(droplet_position - start_position) }));

                if (this->translation_duplication)
                {
                    const auto relative_position = droplet_position - position;

                    this->summary->insert(std::make_shared<geometry::point<float_t>>(position - relative_position),
                        data::polydata_element<float_t, 3>(data_name, data::topology_t::POINT_DATA,
                            { Eigen::Matrix<float_t, 3, 1>(droplet_position - start_position) }));
                }
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_paths(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Create path segments
            const std::string data_name("Time ID");

            this->paths->add(std::make_shared<data::array<std::size_t, 1>>(data_name), data::topology_t::POINT_DATA);

            for (const auto& droplet : all_droplets)
            {
                Eigen::Matrix<float_t, 3, 1> droplet_position = droplet.first.front().position
                    + droplet.first.front().radius * droplet.first.front().translation.normalized();

                for (std::size_t i = 0; i < droplet.first.size(); ++i)
                {
                    const auto next_droplet_position = droplet_position + timesteps[i] * droplet.first[i].translation;

                    this->paths->insert(std::make_shared<geometry::line<float_t>>(
                        geometry::point<float_t>(droplet_position), geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i + 1 }));

                    droplet_position = next_droplet_position;
                }
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_axes(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Create rotation axes
            const std::string data_name("Time ID");

            this->axes->add(std::make_shared<data::array<std::size_t, 1>>(data_name), data::topology_t::CELL_DATA);

            for (const auto& droplet : all_droplets)
            {
                const auto droplet_center = droplet.first.front().position;

                for (std::size_t i = 0; i < droplet.first.size(); ++i)
                {
                    const auto axis_origin = droplet_center - (this->fix_axis_size
                        ? (droplet.first[i].radius * droplet.first[i].rotation.normalized())
                        : (0.5 * this->axis_scale * timesteps[i] * droplet.first[i].rotation));
                    const auto axis_target = droplet_center + (this->fix_axis_size
                        ? (droplet.first[i].radius * droplet.first[i].rotation.normalized())
                        : (0.5 * this->axis_scale * timesteps[i] * droplet.first[i].rotation));

                    this->axes->insert(std::make_shared<geometry::line<float_t>>(
                        geometry::point<float_t>(axis_origin), geometry::point<float_t>(axis_target)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::CELL_DATA, { i }));
                }
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_ribbons(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Create rotation ribbons
            const std::string data_name("Time ID");

            this->ribbons->add(std::make_shared<data::array<std::size_t, 1>>(data_name), data::topology_t::POINT_DATA);

            const auto epsilon_scalar = static_cast<float_t>(1.001);
            const auto epsilon_scalar_outer = epsilon_scalar + this->ribbon_thickness;

            for (const auto& droplet : all_droplets)
            {
                if (droplet.first.front().rotation.isZero()) continue;

                data::polydata<float_t> ribbon;
                ribbon.add(std::make_shared<data::array<std::size_t, 1>>(data_name), data::topology_t::POINT_DATA);

                const auto& initial_droplet_position = droplet.first.front().position;
                const auto& initial_droplet_radius = droplet.first.front().radius;

                Eigen::Matrix<float_t, 3, 1> radial_position = initial_droplet_radius * math::orthonormal(droplet.first.front().rotation).first.normalized();

                Eigen::Matrix<float_t, 3, 1> particle_position = initial_droplet_position + epsilon_scalar * radial_position;
                Eigen::Matrix<float_t, 3, 1> particle_position_outer = initial_droplet_position + epsilon_scalar_outer * radial_position;

                Eigen::Matrix<float_t, 3, 1> pos_min, pos_max, pos_min_outer, pos_max_outer;

                double sum_angle = 0.0;

                for (std::size_t i = 0; i < droplet.first.size(); ++i)
                {
                    math::quaternion<float_t> q;
                    q.from_axis(droplet.first[i].rotation, timesteps[i]);

                    radial_position = q.rotate(radial_position);

                    sum_angle += droplet.first[i].rotation.norm() * timesteps[i];

                    const auto next_particle_position = initial_droplet_position + epsilon_scalar * radial_position;
                    const auto next_particle_position_outer = initial_droplet_position + epsilon_scalar_outer * radial_position;

                    // Create mesh
                    const auto offset_direction = radial_position.cross(next_particle_position - particle_position).normalized();

                    if (i == 0)
                    {
                        pos_min = particle_position - this->ribbon_scale * initial_droplet_radius * offset_direction;
                        pos_max = particle_position + this->ribbon_scale * initial_droplet_radius * offset_direction;

                        pos_min_outer = particle_position_outer - this->ribbon_scale * initial_droplet_radius * offset_direction;
                        pos_max_outer = particle_position_outer + this->ribbon_scale * initial_droplet_radius * offset_direction;

                        // Create mesh for cap
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(pos_min),
                            geometry::point<float_t>(pos_min_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i, i }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(pos_min_outer),
                            geometry::point<float_t>(pos_max_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i, i }));
                    }

                    const auto next_pos_min = next_particle_position - this->ribbon_scale * initial_droplet_radius * offset_direction;
                    const auto next_pos_max = next_particle_position + this->ribbon_scale * initial_droplet_radius * offset_direction;

                    const auto next_pos_min_outer = next_particle_position_outer - this->ribbon_scale * initial_droplet_radius * offset_direction;
                    const auto next_pos_max_outer = next_particle_position_outer + this->ribbon_scale * initial_droplet_radius * offset_direction;

                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(next_pos_max),
                        geometry::point<float_t>(pos_max),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i + 1, i, i }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(next_pos_min),
                        geometry::point<float_t>(next_pos_max),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i + 1, i + 1, i }));

                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_min_outer),
                        geometry::point<float_t>(pos_max_outer),
                        geometry::point<float_t>(next_pos_max_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i, i + 1 }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_min_outer),
                        geometry::point<float_t>(next_pos_max_outer),
                        geometry::point<float_t>(next_pos_min_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i + 1, i + 1 }));

                    {
                        // Create mesh for caps
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(next_pos_max_outer),
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(next_pos_max)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i + 1, i, i + 1 }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_max_outer),
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(next_pos_max_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i, i + 1 }));

                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(next_pos_min),
                            geometry::point<float_t>(pos_min),
                            geometry::point<float_t>(next_pos_min_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i + 1, i, i + 1 }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(next_pos_min_outer),
                            geometry::point<float_t>(pos_min),
                            geometry::point<float_t>(pos_min_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i + 1, i, i }));
                    }

                    particle_position = next_particle_position;
                    particle_position_outer = next_particle_position_outer;

                    pos_min = next_pos_min;
                    pos_max = next_pos_max;

                    pos_min_outer = next_pos_min_outer;
                    pos_max_outer = next_pos_max_outer;
                }

                // Create arrow tip
                const auto tip_rotation = static_cast<float_t>(0.1);
                const auto tip_resolution = 10;
                const auto rotation_increment = tip_rotation / tip_resolution;

                sum_angle += tip_rotation;

                math::quaternion<float_t> q;
                q.from_axis(droplet.first.back().rotation.normalized(), rotation_increment);

                for (std::size_t iteration = 0; iteration < tip_resolution - 1; ++iteration)
                {
                    radial_position = q.rotate(radial_position);

                    const auto next_particle_position = initial_droplet_position + epsilon_scalar * radial_position;
                    const auto next_particle_position_outer = initial_droplet_position + epsilon_scalar_outer * radial_position;

                    const auto offset_direction = radial_position.cross(next_particle_position - particle_position).normalized();

                    const auto scale = (static_cast<float_t>(tip_resolution - iteration) / tip_resolution) * 2 * this->ribbon_scale;
                    const auto next_scale = (static_cast<float_t>(tip_resolution - iteration - 1) / tip_resolution) * 2 * this->ribbon_scale;

                    pos_min = particle_position - scale * initial_droplet_radius * offset_direction;
                    pos_max = particle_position + scale * initial_droplet_radius * offset_direction;

                    pos_min_outer = particle_position_outer - scale * initial_droplet_radius * offset_direction;
                    pos_max_outer = particle_position_outer + scale * initial_droplet_radius * offset_direction;

                    const auto next_pos_min = next_particle_position - next_scale * initial_droplet_radius * offset_direction;
                    const auto next_pos_max = next_particle_position + next_scale * initial_droplet_radius * offset_direction;

                    const auto next_pos_min_outer = next_particle_position_outer - next_scale * initial_droplet_radius * offset_direction;
                    const auto next_pos_max_outer = next_particle_position_outer + next_scale * initial_droplet_radius * offset_direction;

                    if (iteration == 0)
                    {
                        // Create mesh for cap
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_min),
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(pos_max_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_min),
                            geometry::point<float_t>(pos_max_outer),
                            geometry::point<float_t>(pos_min_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    }

                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(next_pos_max),
                        geometry::point<float_t>(pos_max),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(next_pos_min),
                        geometry::point<float_t>(next_pos_max),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_min_outer),
                        geometry::point<float_t>(pos_max_outer),
                        geometry::point<float_t>(next_pos_max_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_min_outer),
                        geometry::point<float_t>(next_pos_max_outer),
                        geometry::point<float_t>(next_pos_min_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                    {
                        // Create mesh for caps
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(next_pos_max),
                            geometry::point<float_t>(next_pos_max_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_max),
                            geometry::point<float_t>(next_pos_max_outer),
                            geometry::point<float_t>(pos_max_outer)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(next_pos_min_outer),
                            geometry::point<float_t>(next_pos_min),
                            geometry::point<float_t>(pos_min)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                        ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(pos_min_outer),
                            geometry::point<float_t>(next_pos_min_outer),
                            geometry::point<float_t>(pos_min)),
                            data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                                { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    }

                    particle_position = next_particle_position;
                    particle_position_outer = next_particle_position_outer;

                    pos_min = next_pos_min;
                    pos_max = next_pos_max;

                    pos_min_outer = next_pos_min_outer;
                    pos_max_outer = next_pos_max_outer;
                }

                q.from_axis(droplet.first.back().rotation.normalized(), rotation_increment);

                radial_position = q.rotate(radial_position);

                const Eigen::Matrix<float_t, 3, 1> next_particle_position = initial_droplet_position + epsilon_scalar * radial_position;
                const Eigen::Matrix<float_t, 3, 1> next_particle_position_outer = initial_droplet_position + epsilon_scalar_outer * radial_position;

                ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                    geometry::point<float_t>(next_particle_position),
                    geometry::point<float_t>(pos_max),
                    geometry::point<float_t>(pos_min)),
                    data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                        { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                    geometry::point<float_t>(pos_min_outer),
                    geometry::point<float_t>(pos_max_outer),
                    geometry::point<float_t>(next_particle_position_outer)),
                    data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                        { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                {
                    // Create mesh for caps
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_max),
                        geometry::point<float_t>(next_particle_position),
                        geometry::point<float_t>(next_particle_position_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_max),
                        geometry::point<float_t>(next_particle_position_outer),
                        geometry::point<float_t>(pos_max_outer)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));

                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(next_particle_position_outer),
                        geometry::point<float_t>(next_particle_position),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                    ribbon.insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(pos_min_outer),
                        geometry::point<float_t>(next_particle_position_outer),
                        geometry::point<float_t>(pos_min)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA,
                            { droplet.first.size(), droplet.first.size(), droplet.first.size() }));
                }

                this->ribbons->merge(ribbon);

                // Duplicate ribbons if they are small enough and the option is set accordingly
                std::vector<float_t> duplication_angles;

                switch (this->ribbon_duplication)
                {
                case dynamic_droplets_aux::ribbon_duplication_t::four:
                    if (sum_angle < (tpf::math::pi<float_t> / 2.0))
                    {
                        duplication_angles.push_back(tpf::math::pi<float_t> / 2.0);
                        duplication_angles.push_back(3.0 * tpf::math::pi<float_t> / 2.0);
                    }

                case dynamic_droplets_aux::ribbon_duplication_t::two:
                    if (sum_angle < tpf::math::pi<float_t>)
                    {
                        duplication_angles.push_back(tpf::math::pi<float_t>);
                    }
                }

                if (!duplication_angles.empty())
                {
                    Eigen::Matrix<float_t, 4, 4> mat_to_origin;
                    mat_to_origin.setIdentity();
                    mat_to_origin.col(3).head(3) -= initial_droplet_position;

                    Eigen::Matrix<float_t, 4, 4> mat_from_origin;
                    mat_from_origin.setIdentity();
                    mat_from_origin.col(3).head(3) += initial_droplet_position;

                    for (const auto& duplication_angle : duplication_angles)
                    {
                        q.from_axis(droplet.first.back().rotation.normalized(), duplication_angle);

                        Eigen::Matrix<float_t, 4, 4> mat_rotation;
                        mat_rotation.setIdentity();
                        mat_rotation.block(0, 0, 3, 3) = q.to_matrix();

                        this->ribbons->merge(*ribbon.clone(tpf::math::transformer<float_t, 3>(mat_from_origin * mat_rotation* mat_to_origin)));
                    }
                }
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_rotation_paths(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Track rotating particle
            const std::string data_name("Time ID");

            this->rotation_paths->add(std::make_shared<data::array<std::size_t, 1>>(data_name), data::topology_t::POINT_DATA);

            for (const auto& droplet : all_droplets)
            {
                if (droplet.first.front().rotation.isZero()) continue;

                Eigen::Matrix<float_t, 3, 1> particle_position = droplet.first.front().position
                    + droplet.first.front().radius * math::orthonormal(droplet.first.front().rotation).first.normalized();

                for (std::size_t i = 0; i < droplet.first.size(); ++i)
                {
                    math::quaternion<float_t> q;
                    q.from_axis(droplet.first[i].rotation, timesteps[i]);

                    Eigen::Matrix<float_t, 3, 1> next_particle_position =
                        q.rotate(Eigen::Matrix<float_t, 3, 1>(particle_position - droplet.first.front().position));

                    next_particle_position = droplet.first.front().position +
                        droplet.first.front().radius * next_particle_position.normalized();

                    this->rotation_paths->insert(std::make_shared<geometry::line<float_t>>(
                        geometry::point<float_t>(particle_position), geometry::point<float_t>(next_particle_position)),
                        data::polydata_element<std::size_t, 1>(data_name, data::topology_t::POINT_DATA, { i, i + 1 }));

                    particle_position = next_particle_position;
                }
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_coordinate_axes(const std::vector<droplet_trace_t>& all_droplets,
            const std::vector<float_t>& timesteps)
        {
            // Transform coordinate system axes
            const std::string data_name("Time information");

            this->coordinate_axes->add(std::make_shared<data::array<float_t, 1>>(data_name), data::topology_t::POINT_DATA);

            for (const auto& droplet : all_droplets)
            {
                Eigen::Matrix<float_t, 3, 1> droplet_position = droplet.first.front().position
                    + droplet.first.front().radius * droplet.first.front().translation.normalized();

                Eigen::Matrix<float_t, 3, 1> x_axis;
                Eigen::Matrix<float_t, 3, 1> y_axis = droplet.first.front().rotation.normalized();
                Eigen::Matrix<float_t, 3, 1> z_axis;

                if (y_axis.isZero()) continue;

                std::tie(x_axis, z_axis) = math::orthonormal(y_axis);

                x_axis *= 0.1f * droplet.first.front().radius;
                y_axis *= 0.1f * droplet.first.front().radius;
                z_axis *= 0.1f * droplet.first.front().radius;

                const auto time_delta = 1 / static_cast<float_t>(droplet.first.size() - 1);

                float_t x_time = 0.0;
                float_t y_time = 2.0;
                float_t z_time = 4.0;

                for (std::size_t i = 0; i < droplet.first.size() - 1; ++i)
                {
                    const auto next_droplet_position = droplet_position + timesteps[i] * droplet.first[i].translation;

                    math::quaternion<float_t> q;
                    q.from_axis(droplet.first[i].rotation, timesteps[i]);

                    const auto next_x_axis = q.rotate(x_axis);
                    const auto next_y_axis = q.rotate(y_axis);
                    const auto next_z_axis = q.rotate(z_axis);

                    const auto next_x_time = x_time + time_delta;
                    const auto next_y_time = y_time + time_delta;
                    const auto next_z_time = z_time + time_delta;

                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position),
                        geometry::point<float_t>(droplet_position + x_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { x_time, x_time, next_x_time }));
                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position + x_axis),
                        geometry::point<float_t>(next_droplet_position + next_x_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { x_time, next_x_time, next_x_time }));

                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position),
                        geometry::point<float_t>(droplet_position + y_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { y_time, y_time, next_y_time }));
                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position + y_axis),
                        geometry::point<float_t>(next_droplet_position + next_y_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { y_time, next_y_time, next_y_time }));

                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position),
                        geometry::point<float_t>(droplet_position + z_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { z_time, z_time, next_z_time }));
                    this->coordinate_axes->insert(std::make_shared<geometry::triangle<float_t>>(
                        geometry::point<float_t>(droplet_position + z_axis),
                        geometry::point<float_t>(next_droplet_position + next_z_axis),
                        geometry::point<float_t>(next_droplet_position)),
                        data::polydata_element<float_t, 1>(data_name, data::topology_t::POINT_DATA, { z_time, next_z_time, next_z_time }));

                    droplet_position += timesteps[i] * droplet.first[i].translation;

                    z_axis = next_z_axis;
                    y_axis = next_y_axis;
                    x_axis = next_x_axis;

                    x_time = next_x_time;
                    y_time = next_y_time;
                    z_time = next_z_time;
                }
            }
        }
    }
}
