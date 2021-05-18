#include "tpf_module_dynamic_droplets.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_line.h"
#include "tpf/geometry/tpf_point.h"
#include "tpf/geometry/tpf_triangle.h"

#include "tpf/log/tpf_log.h"

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
        inline void dynamic_droplets<float_t>::set_algorithm_input(const data::polydata<float_t>& droplets)
        {
            this->droplets = &droplets;
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_output(opt_arg<data::polydata<float_t>> paths,
            opt_arg<data::polydata<float_t>> axes, opt_arg<data::polydata<float_t>> ribbons)
        {
            this->paths = get_or_default<data::polydata<float_t>>(paths);
            this->axes = get_or_default<data::polydata<float_t>>(axes);
            this->ribbons = get_or_default<data::polydata<float_t>>(ribbons);
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::set_algorithm_parameters(const std::size_t num_timesteps,
            const float_t timestep, const float_t axis_scale, const bool axis_translation,
            std::string translation_name, std::string rotation_name, std::string radius_name)
        {
            this->num_timesteps = num_timesteps;
            this->timestep = timestep;
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
            if (this->paths == nullptr && this->axes == nullptr && this->ribbons == nullptr)
            {
                log::info_message(__tpf_info_message("No output selected. Nothing to do."));

                return;
            }

            // Collect droplets for all time steps
            std::vector<data::polydata<float_t>> droplets;
            droplets.reserve(this->num_timesteps);
            droplets.push_back(*this->droplets);

            std::vector<float_t> timesteps;
            timesteps.reserve(this->num_timesteps);
            timesteps.push_back(this->timestep);

            for (std::size_t timestep = 0; timestep < this->num_timesteps; ++timestep)
            {
                // Try to load the next time step
                bool has_available_timestep = true;

                try
                {
                    if (timestep < this->num_timesteps - 1)
                    {
                        // Load next time step's data
                        try
                        {
                            auto next_data = (*this->next_time_frame_callback)();

                            // Get complete data: [Time step delta, droplets]
                            timesteps.push_back(next_data.first);
                            droplets.push_back(next_data.second);
                        }
                        catch (const std::exception&)
                        {
                            has_available_timestep = false;
                        }
                    }

                    if (timestep == this->num_timesteps - 1 || !has_available_timestep)
                    {
                        // Reset data-providing modules
                        this->next_time_frame_callback->reset();
                    }
                }
                catch (...)
                {
                    throw std::runtime_error(__tpf_error_message("Error occured while loading next time frame."));
                }

                log::info_message(__tpf_info_message("Completed time step ", (timestep + 1), " of ", this->num_timesteps, "."));

#ifdef __tpf_debug
                log::info_message(__tpf_info_message("  Number of droplets: ", droplets.back().get_geometry().size()));
                log::info_message(__tpf_info_message("  Time step size: ", timesteps.back()));
#endif

                if (timestep == this->num_timesteps - 1)
                {
                    log::info_message(__tpf_info_message("Completed loading all time steps."));
                }
            }

            // Create output
            if (this->paths != nullptr)
            {
                create_paths(droplets, timesteps);
            }

            if (this->axes != nullptr)
            {
                create_axes(droplets, timesteps);
            }

            if (this->ribbons != nullptr)
            {
                create_ribbons(droplets, timesteps);
            }
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_paths(const std::vector<data::polydata<float_t>>& all_droplets, const std::vector<float_t>& timesteps)
        {
            // Create path segments
            std::size_t num_paths = 0;
            std::vector<std::vector<std::shared_ptr<geometry::line<float_t>>>> paths;

            std::size_t time_index = 0;

            for (const auto& droplets : all_droplets)
            {
                const auto droplet_positions = droplets.get_geometry();
                const auto droplet_translations = *droplets.template get_point_data_as<float_t, 3>(this->translation_name);

                paths.emplace_back(droplet_positions.size());

                for (std::size_t i = 0; i < droplet_positions.size(); ++i)
                {
                    const auto origin = droplet_positions[i]->get_points()[0];
                    const auto target = origin + timesteps[time_index] * droplet_translations(i);

                    paths.back()[i] = std::make_shared<geometry::line<float_t>>(geometry::point<float_t>(origin), geometry::point<float_t>(target));
                }

                num_paths += droplet_positions.size();

                ++time_index;
            }

            // Create array
            auto path_ids_ref = std::make_shared<data::array<std::size_t, 1>>("Time ID", num_paths);
            auto& path_ids = *path_ids_ref;

            // Create map of last targets, which act as new path origins and are initialized
            // as an offset from the droplet center of mass by the radius of the bounding sphere
            std::map<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>> map_to_last_target, next_map_to_last_target;

            const auto droplet_positions = all_droplets.front().get_geometry();
            const auto droplet_translations = *all_droplets.front().template get_point_data_as<float_t, 3>(this->translation_name);
            const auto droplet_radii = *all_droplets.front().template get_point_data_as<float_t, 1>(this->radius_name);

            for (std::size_t i = 0; i < droplet_positions.size(); ++i)
            {
                const auto offset_origin = droplet_positions[i]->get_points()[0] + droplet_radii(i) * droplet_translations(i).normalized();
                const auto path_origin = paths.front()[i]->get_points()[0];

                map_to_last_target[path_origin] = offset_origin;
            }

            // Store path segments (after translation of the segment's origin to the last path segment's target vertex) with IDs
            time_index = 0;
            std::size_t axis_index = 0;

            for (auto it_timestep = paths.cbegin(); it_timestep != paths.cend(); ++it_timestep, ++time_index)
            {
                for (auto it_path = it_timestep->cbegin(); it_path != it_timestep->cend(); ++it_path)
                {
                    const auto points = (*it_path)->get_points();
                    const auto origin = points[0];
                    const auto target = points[1];

                    Eigen::Matrix<float_t, 3, 1> closest = map_to_last_target.cbegin()->first;

                    for (const auto& mapping : map_to_last_target)
                    {
                        if ((closest - origin).norm() > (mapping.first - origin).norm())
                        {
                            closest = mapping.first;
                        }
                    }

                    geometry::point<float_t> new_origin(map_to_last_target[closest]);
                    geometry::point<float_t> new_target(target + (map_to_last_target[closest] - origin));

                    this->paths->insert(std::make_shared<geometry::line<float_t>>(new_origin, new_target));

                    next_map_to_last_target[target] = new_target;

                    path_ids(axis_index++) = time_index;
                }

                std::swap(map_to_last_target, next_map_to_last_target);
                next_map_to_last_target.clear();
            }

            this->paths->add(path_ids_ref, data::topology_t::CELL_DATA);
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_axes(const std::vector<data::polydata<float_t>>& all_droplets, const std::vector<float_t>& timesteps)
        {
            // Create rotation axes
            std::size_t num_axes = 0;
            std::vector<std::vector<std::shared_ptr<geometry::line<float_t>>>> rotation_axes;

            std::size_t time_index = 0;

            for (const auto& droplets : all_droplets)
            {
                const auto droplet_positions = droplets.get_geometry();
                const auto droplet_rotation = *droplets.template get_point_data_as<float_t, 3>(this->rotation_name);

                rotation_axes.emplace_back(droplet_positions.size());

                for (std::size_t i = 0; i < droplet_positions.size(); ++i)
                {
                    const auto center = droplet_positions[i]->get_points()[0];
                    const auto vector = droplet_rotation(i);

                    const auto origin = center - 0.5 * this->axis_scale * timesteps[time_index] * vector;
                    const auto target = center + 0.5 * this->axis_scale * timesteps[time_index] * vector;

                    rotation_axes.back()[i] = std::make_shared<geometry::line<float_t>>(geometry::point<float_t>(origin), geometry::point<float_t>(target));
                }

                num_axes += droplet_positions.size();

                ++time_index;
            }

            // Create array
            auto axes_ids_ref = std::make_shared<data::array<std::size_t, 1>>("Time ID", num_axes);
            auto& axes_ids = *axes_ids_ref;

            // Create map for mapping the axes to their original locations
            std::map<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>> map_to_origin, next_map_to_origin;

            for (auto it_axis = rotation_axes.front().cbegin(); it_axis != rotation_axes.front().cend(); ++it_axis)
            {
                const auto points = (*it_axis)->get_points();
                const auto origin = 0.5 * (points[0] + points[1]);

                map_to_origin[origin] = origin;
            }

            // Store rotation axes (after translation to the droplet origin) with IDs
            time_index = 0;
            std::size_t axis_index = 0;

            for (auto it_timestep = rotation_axes.cbegin(); it_timestep != rotation_axes.cend(); ++it_timestep, ++time_index)
            {
                for (auto it_axis = it_timestep->cbegin(); it_axis != it_timestep->cend(); ++it_axis)
                {
                    if (this->axis_translation)
                    {
                        const auto points = (*it_axis)->get_points();
                        const auto center = 0.5 * (points[0] + points[1]);

                        Eigen::Matrix<float_t, 3, 1> closest = map_to_origin.cbegin()->first;

                        for (const auto& mapping : map_to_origin)
                        {
                            if ((closest - center).norm() > (mapping.first - center).norm())
                            {
                                closest = mapping.first;
                            }
                        }

                        geometry::point<float_t> new_origin(points[0] + (map_to_origin[closest] - center));
                        geometry::point<float_t> new_target(points[1] + (map_to_origin[closest] - center));

                        this->axes->insert(std::make_shared<geometry::line<float_t>>(new_origin, new_target));

                        next_map_to_origin[center] = map_to_origin[closest];
                    }
                    else
                    {
                        this->axes->insert(*it_axis);
                    }

                    axes_ids(axis_index++) = time_index;
                }

                std::swap(map_to_origin, next_map_to_origin);
                next_map_to_origin.clear();
            }

            this->axes->add(axes_ids_ref, data::topology_t::CELL_DATA);
        }

        template <typename float_t>
        inline void dynamic_droplets<float_t>::create_ribbons(const std::vector<data::polydata<float_t>>& all_droplets, const std::vector<float_t>& timesteps)
        {
            // Create rotation axes
            std::vector<std::vector<std::pair<std::shared_ptr<geometry::line<float_t>>, std::shared_ptr<geometry::line<float_t>>>>> rotation_axes;

            std::size_t time_index = 0;

            for (const auto& droplets : all_droplets)
            {
                const auto droplet_positions = droplets.get_geometry();
                const auto droplet_translation = *droplets.template get_point_data_as<float_t, 3>(this->translation_name);
                const auto droplet_rotation = *droplets.template get_point_data_as<float_t, 3>(this->rotation_name);

                rotation_axes.emplace_back(droplet_positions.size());

                for (std::size_t i = 0; i < droplet_positions.size(); ++i)
                {
                    const auto origin = droplet_positions[i]->get_points()[0];
                    const auto target = origin + timesteps[time_index] * droplet_translation(i);
                    const auto vector = droplet_rotation(i);

                    const auto axis_origin = origin - 0.5 * this->axis_scale * timesteps[time_index] * vector;
                    const auto axis_target = origin + 0.5 * this->axis_scale * timesteps[time_index] * vector;

                    rotation_axes.back()[i].first = std::make_shared<geometry::line<float_t>>(geometry::point<float_t>(origin), geometry::point<float_t>(target));
                    rotation_axes.back()[i].second = std::make_shared<geometry::line<float_t>>(geometry::point<float_t>(axis_origin), geometry::point<float_t>(axis_target));
                }

                ++time_index;
            }

            // Create map of last targets, which act as new path origins and are initialized
            // as an offset from the droplet center of mass by the radius of the bounding sphere
            std::map<Eigen::Matrix<float_t, 3, 1>, Eigen::Matrix<float_t, 3, 1>> map_to_last_target, next_map_to_last_target;
            std::map<Eigen::Matrix<float_t, 3, 1>, std::vector<std::pair<std::size_t, geometry::line<float_t>>>>
                map_to_origin_with_corresponding_axes, next_map_to_origin_with_corresponding_axes;

            const auto droplet_positions = all_droplets.front().get_geometry();
            const auto droplet_translations = *all_droplets.front().template get_point_data_as<float_t, 3>(this->translation_name);
            const auto droplet_radii = *all_droplets.front().template get_point_data_as<float_t, 1>(this->radius_name);

            for (std::size_t i = 0; i < droplet_positions.size(); ++i)
            {
                const auto offset_origin = droplet_positions[i]->get_points()[0] + droplet_radii(i) * droplet_translations(i).normalized();
                const auto path_origin = rotation_axes.front()[i].first->get_points()[0];

                map_to_last_target[path_origin] = offset_origin;
                map_to_origin_with_corresponding_axes[path_origin];
            }

            // Store corresponding rotation axes (after translation of the segment's origin to the last path segment's target vertex)
            time_index = 0;
            std::size_t axis_index = 0;

            for (auto it_timestep = rotation_axes.cbegin(); it_timestep != rotation_axes.cend(); ++it_timestep, ++time_index)
            {
                for (auto it_axis = it_timestep->cbegin(); it_axis != it_timestep->cend(); ++it_axis)
                {
                    const auto path = it_axis->first->get_points();
                    const auto axis = it_axis->second->get_points();

                    const auto origin = path[0];
                    const auto target = path[1];

                    const auto axis_origin = axis[0];
                    const auto axis_target = axis[1];
                    const auto vector = axis_target - axis_origin;

                    Eigen::Matrix<float_t, 3, 1> closest = map_to_last_target.cbegin()->first;

                    for (const auto& mapping : map_to_last_target)
                    {
                        if ((closest - origin).norm() > (mapping.first - origin).norm())
                        {
                            closest = mapping.first;
                        }
                    }

                    geometry::point<float_t> new_axis_origin(map_to_last_target[closest] - 0.5 * vector);
                    geometry::point<float_t> new_axis_target(map_to_last_target[closest] + 0.5 * vector);

                    next_map_to_last_target[target] = target + (map_to_last_target[closest] - origin);

                    next_map_to_origin_with_corresponding_axes[target] = map_to_origin_with_corresponding_axes[closest];
                    next_map_to_origin_with_corresponding_axes[target].push_back(std::make_pair(time_index, geometry::line<float_t>(new_axis_origin, new_axis_target)));
                }

                std::swap(map_to_last_target, next_map_to_last_target);
                std::swap(map_to_origin_with_corresponding_axes, next_map_to_origin_with_corresponding_axes);
                next_map_to_last_target.clear();
                next_map_to_origin_with_corresponding_axes.clear();
            }

            // Create array
            const std::size_t num_pieces_per_timestep = 10; // TODO: parameter
            std::size_t num_points = 0;

            for (const auto& axis_map : map_to_origin_with_corresponding_axes)
            {
                num_points += 6 * num_pieces_per_timestep * (axis_map.second.size() - 1);
            }

            auto axes_ids_ref = std::make_shared<data::array<std::size_t, 1>>("Time ID", num_points);
            auto& axes_ids = *axes_ids_ref;

            // Use corresponding axes to generate the ribbon surface
            std::size_t point_index = 0;

            for (const auto& axis_map : map_to_origin_with_corresponding_axes)
            {
                const auto axes = axis_map.second;

                // For all consecutive rotation axes, create triangle mesh
                for (std::size_t axis_index = 0; axis_index < axes.size() - 1; ++axis_index)
                {
                    const auto time_id = axes[axis_index].first;
                    const auto axis = axes[axis_index].second;
                    const auto axis_points = axis.get_points();
                    const auto axis_origin = axis_points[0];
                    const auto axis_target = axis_points[1];
                    const auto axis_direction = axis_target - axis_origin;

                    const auto next_time_id = axes[axis_index + 1].first;
                    const auto next_axis = axes[axis_index + 1].second;
                    const auto next_axis_points = next_axis.get_points();
                    const auto next_axis_origin = next_axis_points[0];
                    const auto next_axis_target = next_axis_points[1];
                    const auto next_axis_direction = next_axis_target - next_axis_origin;

                    const auto step = 1.0f / num_pieces_per_timestep;

                    for (std::size_t step_index = 0; step_index < num_pieces_per_timestep; ++step_index)
                    {
                        const auto current_1 = axis_origin + (step_index * step) * axis_direction;
                        const auto current_2 = axis_origin + ((step_index + 1) * step) * axis_direction;

                        const auto next_1 = next_axis_origin + (step_index * step) * next_axis_direction;
                        const auto next_2 = next_axis_origin + ((step_index + 1) * step) * next_axis_direction;

                        const auto triangle_1 = std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(current_1),
                            geometry::point<float_t>(next_1),
                            geometry::point<float_t>(current_2));

                        const auto triangle_2 = std::make_shared<geometry::triangle<float_t>>(
                            geometry::point<float_t>(current_2),
                            geometry::point<float_t>(next_1),
                            geometry::point<float_t>(next_2));

                        this->ribbons->insert(triangle_1);
                        this->ribbons->insert(triangle_2);

                        axes_ids(point_index++) = time_id;
                        axes_ids(point_index++) = next_time_id;
                        axes_ids(point_index++) = time_id;

                        axes_ids(point_index++) = time_id;
                        axes_ids(point_index++) = next_time_id;
                        axes_ids(point_index++) = next_time_id;
                    }
                }
            }

            this->ribbons->add(axes_ids_ref, data::topology_t::POINT_DATA);
        }
    }
}
