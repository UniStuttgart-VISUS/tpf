#pragma once

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/math/tpf_matrix.h"
#include "tpf/math/tpf_vector.h"

#include <tuple>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace droplets_aux
        {
            /// <summary>
            /// Store cluster information
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class cluster
            {
            public:
                /// <summary>
                /// Constructor
                /// </summary>
                cluster();

                /// <summary>
                /// Add cell to the cluster
                /// </summary>
                /// <param name="coordinates">Coordinates of the cell</param>
                /// <param name="volume">Volume of the cell</param>
                /// <param name="position">Position of the cell</param>
                /// <param name="velocity">Velocity of the cell</param>
                void add(const data::coords3_t& coordinates, floatp_t volume, const math::vec3_t<floatp_t>& position, const math::vec3_t<floatp_t>& velocity);

                /// <summary>
                /// Add cell to the cluster
                /// </summary>
                /// <param name="coordinates">Coordinates of the cell</param>
                /// <param name="volume">Volume of the cell</param>
                /// <param name="position">Position of the cell</param>
                void add(const data::coords3_t& coordinates, floatp_t volume, const math::vec3_t<floatp_t>& position);

                /// <summary>
                /// Add the contents of another cluster to this one
                /// </summary>
                /// <param name="other">Other cluster to copy from</param>
                /// <returns>This</returns>
                cluster<floatp_t>& operator+=(const cluster<floatp_t>& other);

                /// <summary>
                /// Return the coordinates associated with this cluster
                /// </summary>
                /// <returns>Associated coordinates</returns>
                const std::vector<data::coords3_t>& get_coordinates() const;

                /// <summary>
                /// Return the position (barycenter) of this cluster
                /// </summary>
                /// <returns>Barycenter</returns>
                math::vec3_t<floatp_t> get_position() const;

                /// <summary>
                /// Return the average velocity of this cluster
                /// </summary>
                /// <returns>Average velocity and its energy</returns>
                std::pair<math::vec3_t<floatp_t>, floatp_t> get_velocity() const;

                /// <summary>
                /// Return the volume this cluster
                /// </summary>
                /// <returns>Volume</returns>
                floatp_t get_volume() const;

                /// <summary>
                /// Return the kinetic energy of this cluster
                /// </summary>
                /// <returns>Kinetic energy</returns>
                floatp_t get_energy() const;

                /// <summary>
                /// Return the rotation of this cluster, using inertia and angular momentum
                /// </summary>
                /// <returns>Rotation axis, its rotation energy and resulting local energy</returns>
                std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> get_rotation_mechanics() const;

                /// <summary>
                /// Return the rotation of this cluster, minimizing velocities
                /// </summary>
                /// <returns>Rotation axis, its rotation energy and resulting local energy</returns>
                std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> get_rotation_velocities() const;

                /// <summary>
                /// Return the rotation of this cluster, using PCA
                /// </summary>
                /// <returns>Rotation axis, its rotation energy and resulting local energy</returns>
                std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> get_rotation_pca() const;

                /// <summary>
                /// Return the inertia of this cluster
                /// </summary>
                /// <returns>Inertia</returns>
                math::mat3_t<floatp_t> get_inertia() const;

                /// <summary>
                /// Return the angular momentum of this cluster
                /// </summary>
                /// <returns>Angular momentum</returns>
                math::vec3_t<floatp_t> get_angular_momentum() const;

                /// <summary>
                /// Return a sphere radius that completely envelops the droplet
                /// </summary>
                /// <returns>Bounding sphere radius</returns>
                floatp_t get_bounding_sphere() const;

            private:
                /// <summary>
                /// Calculate the rotational and local energy
                /// </summary>
                /// <returns>Rotation energy and resulting local energy</returns>
                std::pair<floatp_t, floatp_t> calculate_energies(const math::vec3_t<floatp_t>& rotation_axis) const;

                /// Coordinates of cells in this cluster
                std::vector<data::coords3_t> coordinates;

                /// Positions and velocities of cells in this cluster
                std::vector<math::vec3_t<floatp_t>> positions;
                std::vector<math::vec3_t<floatp_t>> velocities;

                /// Sum of the positions and the velocities
                math::vec3_t<floatp_t> sum_position;
                math::vec3_t<floatp_t> sum_velocity;

                /// Volumes of cells in this cluster
                std::vector<floatp_t> volumes;

                /// Volume and total energy of the cluster
                floatp_t volume;
                floatp_t energy;
            };
        }
    }
}

#include "tpf_cluster.inl"
