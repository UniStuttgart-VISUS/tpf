#include "tpf_module_droplets.h"

#include "tpf_cluster.h"

#include "tpf/data/tpf_array.h"
#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"
#include "tpf/data/tpf_polydata.h"

#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_linalg.h"
#include "tpf/math/tpf_vector.h"

#include "tpf/mpi/tpf_mpi_grid.h"

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline std::size_t droplets<float_t>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t>
        inline droplets<float_t>::droplets() { }

        template <typename float_t>
        inline std::string droplets<float_t>::get_name() const
        {
            return std::string("Droplets");
        }

        template <typename float_t>
        inline void droplets<float_t>::set_algorithm_input(const data::grid<float_t, float_t, 3, 1>& fractions, const data::grid<float_t, float_t, 3, 3>& positions,
            std::optional<std::reference_wrapper<const data::grid<float_t, float_t, 3, 3>>> velocities)
        {
            this->fractions = &fractions;
            this->positions = &positions;
            this->velocities = get_or_default<const data::grid<float_t, float_t, 3, 3>>(velocities);
        }

        template <typename float_t>
        inline void droplets<float_t>::set_algorithm_output(data::grid<long long, float_t, 3, 1>& droplet_ids, data::grid<float_t, float_t, 3, 1>& droplet_volumes,
            std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>> droplet_velocities, data::polydata<float_t>& droplets)
        {
            this->droplet_ids = &droplet_ids;
            this->droplet_volumes = &droplet_volumes;
            this->droplet_velocities = get_or_default<data::grid<float_t, float_t, 3, 3>>(droplet_velocities);
            this->local_droplets = &droplets;
        }

        template <typename float_t>
        inline void droplets<float_t>::set_algorithm_parameters(bool calculate_translation, bool calculate_rotation, bool calculate_energy, bool calculate_inertia,
            droplets_aux::scale_method_t scaling_method, droplets_aux::rotation_method_t rotation_method)
        {
            this->calculate_translation = calculate_translation;
            this->calculate_rotation = calculate_rotation;
            this->calculate_energy = calculate_energy;
            this->calculate_inertia = calculate_inertia;

            this->scale_method = scaling_method;
            this->rotation_method = rotation_method;
        }

        template <typename float_t>
        inline void droplets<float_t>::run_algorithm()
        {
            // Sanity checks
            if (this->velocities == nullptr)
            {
                if (this->calculate_translation || this->calculate_rotation || this->calculate_energy || this->calculate_inertia)
                {
                    log::warning_message(__tpf_warning_message("No velocity field given. Cannot calculate translation, rotation, energy and inertia."));
                }

                this->calculate_translation = this->calculate_rotation = this->calculate_energy = this->calculate_inertia = false;
            }
            else if (!this->calculate_translation && this->calculate_rotation)
            {
                log::warning_message(__tpf_warning_message("Calculation of rotation without prior calculation of translation not allowed. Calculating both instead."));

                this->calculate_translation = this->calculate_rotation = true;
            }

            // Get global scope of input
            const auto complete_fractions = mpi::mpi_grid::all_gather(*this->fractions);
            const auto complete_positions = mpi::mpi_grid::all_gather(*this->positions);
            const auto complete_velocities = mpi::mpi_grid::all_gather((this->velocities != nullptr) ? *this->velocities : data::grid<float_t, float_t, 3, 3>());

            // Create and initialize complete droplet IDs and droplets
            data::grid<long long, float_t, 3, 1> complete_droplet_indices("Droplet IDs", complete_fractions.get_extent());
            complete_droplet_indices.initialize(-1LL);

            data::polydata<float_t> complete_droplets;

            // Clustering on neighborhood, generating droplets
            std::map<long long, droplets_aux::cluster<float_t>> complete_clusters;

            long long next_cluster_id = 0LL;

            for (std::size_t z = complete_fractions.get_extent()[2].first; z <= complete_fractions.get_extent()[2].second; ++z)
            {
                for (std::size_t y = complete_fractions.get_extent()[1].first; y <= complete_fractions.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = complete_fractions.get_extent()[0].first; x <= complete_fractions.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        if (complete_fractions(coords) > static_cast<float_t>(0.0L))
                        {
                            // For all neighbors
                            std::vector<long long> clusters_to_merge;

                            for (long long k = -1; k <= 1; ++k)
                            {
                                for (long long j = -1; j <= 1; ++j)
                                {
                                    for (long long i = -1; i <= 1; ++i)
                                    {
                                        const data::coords3_t neighbor_coords = coords + data::coords3_t(i, j, k);

                                        if ((i != 0 || j != 0 || k != 0) && complete_fractions.is_on_grid(neighbor_coords)
                                            && complete_fractions(neighbor_coords) > static_cast<float_t>(0.0L))
                                        {
                                            // If neighbor belongs to a droplet, add this point to it
                                            if (complete_droplet_indices(coords) == -1 && complete_droplet_indices(neighbor_coords) != -1)
                                            {
                                                complete_droplet_indices(coords) = complete_droplet_indices(neighbor_coords);
                                            }

                                            // If there is a different droplet, it needs to be merged into this one
                                            if (complete_droplet_indices(neighbor_coords) != -1)
                                            {
                                                clusters_to_merge.push_back(complete_droplet_indices(neighbor_coords));
                                            }
                                        }
                                    }
                                }
                            }

                            // If there was no droplet, create new one
                            if (clusters_to_merge.size() == 0)
                            {
                                complete_droplet_indices(coords) = next_cluster_id++;
                                complete_clusters[complete_droplet_indices(coords)];
                            }

                            // Add influence of current cell
                            const float_t volume_or_mass = complete_fractions(coords) * complete_fractions.get_cell_sizes(coords).prod();

                            if (this->velocities != nullptr)
                            {
                                complete_clusters.at(complete_droplet_indices(coords)).add(coords, volume_or_mass, complete_positions(coords), complete_velocities(coords));
                            }
                            else
                            {
                                complete_clusters.at(complete_droplet_indices(coords)).add(coords, volume_or_mass, complete_positions(coords));
                            }

                            // If there are multiple droplets, merge them
                            if (clusters_to_merge.size() >= 2)
                            {
                                for (std::size_t i = 1; i < clusters_to_merge.size(); ++i)
                                {
                                    if (clusters_to_merge[0] != clusters_to_merge[i] && complete_clusters.find(clusters_to_merge[i]) != complete_clusters.end())
                                    {
                                        for (const auto& coord : complete_clusters.at(clusters_to_merge[i]).get_coordinates())
                                        {
                                            complete_droplet_indices(coord) = clusters_to_merge[0];
                                        }

                                        complete_clusters.at(clusters_to_merge[0]) += complete_clusters.at(clusters_to_merge[i]);

                                        complete_clusters.erase(clusters_to_merge[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Change keys in order to shrink the range of droplets IDs and create a map for the droplet ID grid
            std::map<long long, long long> map_ids;

            next_cluster_id = 0LL;

            for (auto it = complete_clusters.begin(); it != complete_clusters.end(); ++next_cluster_id)
            {
                if (next_cluster_id != it->first)
                {
                    map_ids[it->first] = next_cluster_id;

                    std::swap(complete_clusters[next_cluster_id], it->second);
                    complete_clusters.erase(it++);
                }
                else
                {
                    map_ids[next_cluster_id] = next_cluster_id;
                    ++it;
                }
            }

            // Create an object for each droplet
            for (const auto& cluster : complete_clusters)
            {
                const droplets_aux::cluster<float_t>& cluster_info = cluster.second;

                auto point = std::make_shared<geometry::point<float_t>>(cluster_info.get_position());
                complete_droplets.insert(point);

#ifdef __tpf_debug
                log::info_message(__tpf_info_message("Added droplet with barycenter (", point->get_vertex()[0], ", ", point->get_vertex()[1], ", ", point->get_vertex()[2], ")"));
#endif
            }

            // Create arrays for the poly data
            data::array<long long, 1>& ids = *complete_droplets.template create<long long, 1>("Droplet IDs", data::topology_t::POINT_DATA);
            data::array<float_t, 1>& volume = *complete_droplets.template create<float_t, 1>("Volumes", data::topology_t::POINT_DATA);
            data::array<float_t, 1>& radius = *complete_droplets.template create<float_t, 1>("Radius", data::topology_t::POINT_DATA);

            std::shared_ptr<data::array<float_t, 3>> translation_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 3>> rotation_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 1>> energy_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 1>> energy_translation_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 1>> energy_rotation_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 1>> energy_local_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 3, 3>> inertia_ptr = nullptr;

            if (this->calculate_translation) translation_ptr = complete_droplets.template create<float_t, 3>("Translation", data::topology_t::POINT_DATA);
            if (this->calculate_rotation) rotation_ptr = complete_droplets.template create<float_t, 3>("Rotation", data::topology_t::POINT_DATA);
            if (this->calculate_energy) energy_ptr = complete_droplets.template create<float_t, 1>("Energy", data::topology_t::POINT_DATA);
            if (this->calculate_energy && this->calculate_translation) energy_translation_ptr = complete_droplets.template create<float_t, 1>("Translation Energy", data::topology_t::POINT_DATA);
            if (this->calculate_energy && this->calculate_rotation) energy_rotation_ptr = complete_droplets.template create<float_t, 1>("Rotation Energy", data::topology_t::POINT_DATA);
            if (this->calculate_energy && this->calculate_rotation) energy_local_ptr = complete_droplets.template create<float_t, 1>("Local Energy", data::topology_t::POINT_DATA);
            if (this->calculate_inertia) inertia_ptr = complete_droplets.template create<float_t, 3, 3>("Inertia", data::topology_t::POINT_DATA);

#ifdef __tpf_detailed
            std::shared_ptr<data::array<float_t, 3>> inertia_vec_1_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 3>> inertia_vec_2_ptr = nullptr;
            std::shared_ptr<data::array<float_t, 3>> inertia_vec_3_ptr = nullptr;

            if (this->calculate_inertia)
            {
                inertia_vec_1_ptr = complete_droplets.template create<float_t, 3>("Inertia #1", data::topology_t::POINT_DATA);
                inertia_vec_2_ptr = complete_droplets.template create<float_t, 3>("Inertia #2", data::topology_t::POINT_DATA);
                inertia_vec_3_ptr = complete_droplets.template create<float_t, 3>("Inertia #3", data::topology_t::POINT_DATA);
            }
#endif

            // Droplet-wise calculate different quantities
            for (const auto& cluster : complete_clusters)
            {
                const long long cluster_id = cluster.first;
                const droplets_aux::cluster<float_t>& cluster_info = cluster.second;

                // Store ID, volume and bounding sphere radius
                ids(cluster_id) = cluster_id;
                volume(cluster_id) = cluster_info.get_volume();
                radius(cluster_id) = cluster_info.get_bounding_sphere();

#ifdef __tpf_debug
                log::info_message(__tpf_info_message("Droplet ID:\t", cluster_id));
                log::info_message(__tpf_info_message("Volume:\t", volume(cluster_id)));
                log::info_message(__tpf_info_message("Radius:\t", radius(cluster_id)));
#endif

                // Store translation
                if (this->calculate_translation)
                {
                    const auto velocity = cluster_info.get_velocity();

                    auto& translation = *translation_ptr;
                    translation(cluster_id) = velocity.first;

#ifdef __tpf_debug
                    log::info_message(__tpf_info_message("Translation:\t(", translation(cluster_id)[0], ", ", translation(cluster_id)[1], ", ", translation(cluster_id)[2], ")"));
#endif

                    if (this->calculate_energy)
                    {
                        auto& energy = *energy_translation_ptr;
                        energy(cluster_id) = velocity.second;

#ifdef __tpf_debug
                        log::info_message(__tpf_info_message("Energy (translation):\t", energy(cluster_id)));
#endif
                    }
                }

                // Calculate and store rotation
                if (this->calculate_rotation)
                {
                    std::tuple<math::vec3_t<float_t>, float_t, float_t> rotation_axis;

                    auto& rotation = *rotation_ptr;

                    switch (this->rotation_method)
                    {
                    case droplets_aux::rotation_method_t::MECHANICS:
                        rotation_axis = cluster_info.get_rotation_mechanics();
                        break;
                    case droplets_aux::rotation_method_t::VELOCITIES:
                        rotation_axis = cluster_info.get_rotation_velocities();
                        break;
                    case droplets_aux::rotation_method_t::PCA:
                        rotation_axis = cluster_info.get_rotation_pca();
                        break;
                    }

                    rotation(cluster_id) = std::get<0>(rotation_axis);

#ifdef __tpf_debug
                    log::info_message(__tpf_info_message("Rotation:\t(", rotation(cluster_id)[0], ", ", rotation(cluster_id)[1], ", ", rotation(cluster_id)[2], ")"));
#endif

                    if (this->calculate_energy)
                    {
                        auto& rotation_energy = *energy_rotation_ptr;
                        rotation_energy(cluster_id) = std::get<1>(rotation_axis);

#ifdef __tpf_debug
                        log::info_message(__tpf_info_message("Energy (rotation):\t", rotation_energy(cluster_id)));
#endif

                        auto& local_energy = *energy_local_ptr;
                        local_energy(cluster_id) = std::get<2>(rotation_axis);

#ifdef __tpf_debug
                        log::info_message(__tpf_info_message("Energy (local):\t", local_energy(cluster_id)));
#endif
                    }
                }

                // Store energy
                if (this->calculate_energy)
                {
                    auto& energy = *energy_ptr;
                    energy(cluster_id) = cluster_info.get_energy();

#ifdef __tpf_debug
                    log::info_message(__tpf_info_message("Energy:\t", energy(cluster_id)));
#endif
                }

                // Calculate and store inertia
                if (this->calculate_inertia)
                {
                    auto& inertia = *inertia_ptr;
                    inertia(cluster_id) = cluster_info.get_inertia();

#ifdef __tpf_detailed
                    auto& inertia_vec_1 = *inertia_vec_1_ptr;
                    auto& inertia_vec_2 = *inertia_vec_2_ptr;
                    auto& inertia_vec_3 = *inertia_vec_3_ptr;

                    // Calculate and save Eigenvalues and -vectors
                    const math::eigenpairs<float_t, 3> eigenpairs = math::calculate_eigenpair<float_t, 3>(inertia(cluster_id));

                    // Scale eigenvalues
                    const float_t eigenvalue_1 = scale(eigenpairs.eigenvalues[0]);
                    const float_t eigenvalue_2 = scale(eigenpairs.eigenvalues[1]);
                    const float_t eigenvalue_3 = scale(eigenpairs.eigenvalues[2]);

                    const math::vec3_t<float_t> eigenvector_1 = eigenpairs.eigenvectors[0];
                    const math::vec3_t<float_t> eigenvector_2 = eigenpairs.eigenvectors[1];
                    const math::vec3_t<float_t> eigenvector_3 = eigenpairs.eigenvectors[2];

                    // Write scaled eigenvectors
                    inertia_vec_1(cluster_id) = eigenvector_1 * eigenvalue_1;
                    inertia_vec_2(cluster_id) = eigenvector_2 * eigenvalue_2;
                    inertia_vec_3(cluster_id) = eigenvector_3 * eigenvalue_3;
#endif

#ifdef __tpf_debug
                    log::info_message(__tpf_info_message("Inertia tensor:\t",
                        "\t\t\t((", inertia(cluster_id)(0, 0), ", ", inertia(cluster_id)(0, 1), ", ", inertia(cluster_id)(0, 2), ")\n",
                        "\t\t\t (", inertia(cluster_id)(1, 0), ", ", inertia(cluster_id)(1, 1), ", ", inertia(cluster_id)(1, 2), ")\n",
                        "\t\t\t (", inertia(cluster_id)(2, 0), ", ", inertia(cluster_id)(2, 1), ", ", inertia(cluster_id)(2, 2), "))"));
#endif
                }
            }

            // Store quantities in grids
            const data::grid<float_t, float_t, 3, 1>& fractions = *this->fractions;

            data::grid<long long, float_t, 3, 1>& droplet_ids = *this->droplet_ids;
            data::grid<float_t, float_t, 3, 1>& droplet_volumes = *this->droplet_volumes;
            data::grid<float_t, float_t, 3, 3>& droplet_velocities = *this->droplet_velocities;

            for (std::size_t z = droplet_ids.get_extent()[2].first; z <= droplet_ids.get_extent()[2].second; ++z)
            {
                for (std::size_t y = droplet_ids.get_extent()[1].first; y <= droplet_ids.get_extent()[1].second; ++y)
                {
                    for (std::size_t x = droplet_ids.get_extent()[0].first; x <= droplet_ids.get_extent()[0].second; ++x)
                    {
                        const data::coords3_t coords(x, y, z);

                        // Set initial values
                        droplet_ids(coords) = -1LL;
                        droplet_volumes(coords) = static_cast<float_t>(0.0L);
                        droplet_velocities(coords).setZero();

                        // Store calculated values
                        if (fractions(coords) > static_cast<float_t>(0.0L) && complete_droplet_indices(coords) != -1LL)
                        {
                            const long long cluster_id = map_ids.at(complete_droplet_indices(coords));

                            droplet_ids(coords) = cluster_id;
                            droplet_volumes(coords) = volume(cluster_id);

                            if (this->calculate_translation)
                            {
                                auto& translation = *translation_ptr;
                                droplet_velocities(coords) += translation(cluster_id);
                            }

                            if (this->calculate_rotation)
                            {
                                auto& rotation = *rotation_ptr;
                                droplet_velocities(coords) += rotation(cluster_id).cross(complete_positions(coords) - complete_droplets.get_object(cluster_id)->get_points()[0]);
                            }
                        }
                    }
                }
            }

            // Extract local poly data
            if (mpi::get_instance().check_mpi_status())
            {
                const data::extent_t extent = mpi::mpi_grid::get_local_extent(fractions.get_extent()).first;
                const data::coords3_t min_extent(extent[0].first, extent[1].first, extent[2].first);
                const data::coords3_t max_extent(extent[0].second + 1, extent[1].second + 1, extent[2].second + 1);

                const auto min_coordinates = fractions.get_node_coordinates(min_extent);
                const auto max_coordinates = fractions.get_node_coordinates(max_extent);

                data::area_t<float_t> area(3);
                area[0].first = min_coordinates[0];
                area[1].first = min_coordinates[1];
                area[2].first = min_coordinates[2];
                area[0].second = max_coordinates[0];
                area[1].second = max_coordinates[1];
                area[2].second = max_coordinates[2];

                *this->local_droplets = std::move(complete_droplets.extract_area(area, true));
            }
            else
            {
                *this->local_droplets = std::move(complete_droplets);
            }
        }

        template <typename float_t>
        inline float_t droplets<float_t>::scale(float_t value) const
        {
            switch (this->scale_method)
            {
            case droplets_aux::scale_method_t::NORMALIZED:
                value = static_cast<float_t>(1.0L);

                break;
            case droplets_aux::scale_method_t::BY_LOGARITHMIC_EIGENVALUE:
                value = std::log2(value);
            }

            return value;
        }
    }
}
