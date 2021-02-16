#include "tpf_cluster.h"

#include "tpf/algorithm/tpf_least_squares.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_linalg.h"
#include "tpf/math/tpf_matrix.h"
#include "tpf/math/tpf_vector.h"

#include "Eigen/Dense"

#include <algorithm>
#include <new>
#include <tuple>
#include <utility>
#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace droplets_aux
        {
            template <typename floatp_t>
            inline cluster<floatp_t>::cluster() : volume(0.0), energy(0.0), sum_position(0.0, 0.0, 0.0), sum_velocity(0.0, 0.0, 0.0)
            { }

            template <typename floatp_t>
            inline void cluster<floatp_t>::add(const data::coords3_t& coordinates, const floatp_t volume, const math::vec3_t<floatp_t>& position, const math::vec3_t<floatp_t>& velocity)
            {
                this->coordinates.push_back(coordinates);

                this->positions.push_back(position);
                this->velocities.push_back(velocity);

                this->sum_position += volume * position;
                this->sum_velocity += volume * velocity;

                this->volumes.push_back(volume);

                this->volume += volume;
                this->energy += volume * velocity.squaredNorm();
            }

            template <typename floatp_t>
            inline void cluster<floatp_t>::add(const data::coords3_t& coordinates, const floatp_t volume, const math::vec3_t<floatp_t>& position)
            {
                this->coordinates.push_back(coordinates);

                this->positions.push_back(position);

                this->sum_position += volume * position;

                this->volumes.push_back(volume);

                this->volume += volume;
            }

            template <typename floatp_t>
            inline cluster<floatp_t>& cluster<floatp_t>::operator+=(const cluster<floatp_t>& other)
            {
                this->coordinates.insert(this->coordinates.end(), other.coordinates.begin(), other.coordinates.end());

                this->positions.insert(this->positions.end(), other.positions.begin(), other.positions.end());
                this->velocities.insert(this->velocities.end(), other.velocities.begin(), other.velocities.end());

                this->sum_position += other.sum_position;
                this->sum_velocity += other.sum_velocity;

                this->volumes.insert(this->volumes.end(), other.volumes.begin(), other.volumes.end());

                this->volume += other.volume;
                this->energy += other.energy;

                return *this;
            }

            template <typename floatp_t>
            inline const std::vector<data::coords3_t>& cluster<floatp_t>::get_coordinates() const
            {
                return this->coordinates;
            }

            template <typename floatp_t>
            inline math::vec3_t<floatp_t> cluster<floatp_t>::get_position() const
            {
                return this->sum_position / this->volume;
            }

            template <typename floatp_t>
            inline std::pair<math::vec3_t<floatp_t>, floatp_t> cluster<floatp_t>::get_velocity() const
            {
                return std::make_pair<math::vec3_t<floatp_t>, floatp_t>(this->sum_velocity / this->volume,
                    static_cast<floatp_t>(0.5L) * this->sum_velocity.squaredNorm());
            }

            template <typename floatp_t>
            inline floatp_t cluster<floatp_t>::get_volume() const
            {
                return this->volume;
            }

            template <typename floatp_t>
            inline floatp_t cluster<floatp_t>::get_energy() const
            {
                return static_cast<floatp_t>(0.5L) * this->energy;
            }

            template <typename floatp_t>
            inline std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> cluster<floatp_t>::get_rotation_mechanics() const
            {
                math::vec3_t<floatp_t> rotation_axis;
                rotation_axis.setZero();

                floatp_t rotation_energy = static_cast<floatp_t>(0.0L);
                floatp_t local_energy = static_cast<floatp_t>(0.0L);

                // Compute inertia and angular momentum and use it to solve for the rotation axis
                const std::size_t num_entries = this->positions.size();

                if (num_entries > 2)
                {
                    try
                    {
                        // Compute rotation axis
                        const auto inertia = get_inertia();
                        const auto angular_momentum = get_angular_momentum();

                        rotation_axis = inertia.inverse() * angular_momentum;

                        // Calculate rotational energy and resulting local energy
                        const auto energies = calculate_energies(rotation_axis);

                        rotation_energy = energies.first;
                        local_energy = energies.second;
                    }
                    catch (const std::bad_alloc & ex)
                    {
                        log::warning_message(__tpf_nested_warning_message(ex.what(), "Rotation for large droplet of ", num_entries, " cells could not be calculated."));

                        rotation_axis.setZero();
                    }
                }

                return std::make_tuple(rotation_axis, static_cast<floatp_t>(0.5L)* rotation_energy, static_cast<floatp_t>(0.5L)* local_energy);
            }

            template <typename floatp_t>
            inline std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> cluster<floatp_t>::get_rotation_velocities() const
            {
                math::vec3_t<floatp_t> rotation_axis;
                rotation_axis.setZero();

                floatp_t rotation_energy = static_cast<floatp_t>(0.0L);
                floatp_t local_energy = static_cast<floatp_t>(0.0L);

                const auto barycenter = get_position();
                const auto translation = get_velocity().first;

                // Fill matrices for the least squares problem
                const std::size_t num_entries = this->positions.size();

                if (num_entries > 2)
                {
                    try
                    {
                        Eigen::Matrix<floatp_t, Eigen::Dynamic, 1> y;
                        y.resize(num_entries * 3, Eigen::NoChange);

                        Eigen::Matrix<floatp_t, Eigen::Dynamic, 3> X;
                        X.resize(num_entries * 3, Eigen::NoChange);

                        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
                        {
                            const auto position = this->positions[entry_index] - barycenter;
                            const auto velocity = this->velocities[entry_index] - translation;

                            y.block(entry_index * 3, 0, 3, 1) = velocity;
                            X.block(entry_index * 3, 0, 3, 3) << 0.0, position[2], -position[1], -position[2], 0.0, position[0], position[1], -position[0], 0.0;
                        }

                        // Solve minimization problem
                        rotation_axis = algorithm::template least_squares(X, y);

                        // Calculate rotational energy and resulting local energy
                        const auto energies = calculate_energies(rotation_axis);

                        rotation_energy = energies.first;
                        local_energy = energies.second;
                    }
                    catch (const std::bad_alloc& ex)
                    {
                        log::warning_message(__tpf_nested_warning_message(ex.what(), "Rotation for large droplet of ", num_entries, " cells could not be calculated."));

                        rotation_axis.setZero();
                    }
                }

                return std::make_tuple(rotation_axis, static_cast<floatp_t>(0.5L) * rotation_energy, static_cast<floatp_t>(0.5L) * local_energy);
            }

            template <typename floatp_t>
            inline std::tuple<math::vec3_t<floatp_t>, floatp_t, floatp_t> cluster<floatp_t>::get_rotation_pca() const
            {
                math::vec3_t<floatp_t> rotation_axis;
                rotation_axis.setZero();

                floatp_t rotation_energy = static_cast<floatp_t>(0.0L);
                floatp_t local_energy = static_cast<floatp_t>(0.0L);

                const auto barycenter = get_position();
                const auto translation = get_velocity().first;

                // Fill matrix for the computation of the covariance matrix
                const std::size_t num_entries = this->positions.size();

                if (num_entries > 2)
                {
                    try
                    {
                        Eigen::Matrix<floatp_t, 3, Eigen::Dynamic> X;
                        X.resize(Eigen::NoChange, num_entries);

                        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
                        {
                            const auto position = this->positions[entry_index] - barycenter;

                            X.block(0, entry_index, 3, 1) = position;
                        }

                        // Compute covariance matrix, and eigenvectors and -values
                        const auto C = (1.0 / (num_entries - 1)) * X * X.transpose();

                        const auto eigenpairs = math::calculate_eigenpair<floatp_t, 3>(C);

                        // For each eigenvector, compute as if it was the axis of rotation
                        auto eigenvectors = eigenpairs.eigenvectors;

                        for (auto& axis : eigenvectors)
                        {
                            // Minimize velocities, scaling the axis accordingly
                            std::vector<floatp_t> scalars(num_entries * 3);

                            auto median = [](std::vector<floatp_t>& v) -> floatp_t
                            {
                                const std::size_t n = v.size() / 2;
                                std::nth_element(v.begin(), v.begin() + n, v.end());
                                return v[n];
                            };

                            for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
                            {
                                const auto position = this->positions[entry_index] - barycenter;
                                const auto velocity = this->velocities[entry_index] - translation;

                                const auto values = velocity.cwiseQuotient(axis.cross(position));

                                scalars[entry_index * 3 + 0] = values[0];
                                scalars[entry_index * 3 + 1] = values[1];
                                scalars[entry_index * 3 + 2] = values[2];
                            }

                            axis = axis.normalized() * median(scalars);
                        }

                        // Pick the largest axis
                        rotation_axis = (eigenvectors[0].norm() > eigenvectors[1].norm() && eigenvectors[0].norm() > eigenvectors[2].norm()) ? eigenvectors[0] :
                            ((eigenvectors[1].norm() > eigenvectors[2].norm()) ? eigenvectors[1] : eigenvectors[2]);

                        // Calculate rotational energy
                        const auto energies = calculate_energies(rotation_axis);

                        rotation_energy = energies.first;
                        local_energy = energies.second;
                    }
                    catch (const std::bad_alloc & ex)
                    {
                        log::warning_message(__tpf_nested_warning_message(ex.what(), "Rotation for large droplet of ", num_entries, " cells could not be calculated."));

                        rotation_axis.setZero();
                    }
                }

                return std::make_tuple(rotation_axis, static_cast<floatp_t>(0.5L) * rotation_energy, static_cast<floatp_t>(0.5L) * local_energy);
            }

            template <typename floatp_t>
            inline math::mat3_t<floatp_t> cluster<floatp_t>::get_inertia() const
            {
                Eigen::Matrix<floatp_t, 3, 3> inertia;
                inertia.setZero();

                const auto barycenter = get_position();

                for (std::size_t i = 0; i < this->positions.size(); ++i)
                {
                    const auto local_position = this->positions[i] - barycenter;

                    Eigen::Matrix<floatp_t, 3, 3> local;
                    local.setZero();
                    local.diagonal().setConstant(local_position.squaredNorm());
                    local -= local_position * local_position.transpose();

                    inertia += this->volumes[i] * local;
                }

                return inertia;
            }

            template <typename floatp_t>
            inline math::vec3_t<floatp_t> cluster<floatp_t>::get_angular_momentum() const
            {
                Eigen::Matrix<floatp_t, 3, 1> angular_momentum;
                angular_momentum.setZero();

                const auto barycenter = get_position();
                const auto translation = get_velocity().first;

                for (std::size_t i = 0; i < this->positions.size(); ++i)
                {
                    const auto local_position = this->positions[i] - barycenter;
                    const auto local_velocity = this->velocities[i] - translation;
                    const auto mass = this->volumes[i];

                    angular_momentum += local_position.cross(mass * local_velocity);
                }

                return angular_momentum;
            }

            template <typename floatp_t>
            inline floatp_t cluster<floatp_t>::get_bounding_sphere() const
            {
                // Compute the maximum distance of all positions to the center of mass
                floatp_t radius = 0.0f;

                const auto barycenter = get_position();

                for (std::size_t i = 0; i < this->positions.size(); ++i)
                {
                    const auto local_position = this->positions[i] - barycenter;

                    radius = std::max(radius, local_position.norm());
                }

                // Return bounding sphere radius, which is slightly further away than the remotest point
                return radius * 1.1f;
            }

            template <typename floatp_t>
            std::pair<floatp_t, floatp_t> cluster<floatp_t>::calculate_energies(const math::vec3_t<floatp_t>& rotation_axis) const
            {
                const std::size_t num_entries = this->positions.size();

                const auto barycenter = get_position();
                const auto translation = get_velocity().first;

                floatp_t rotation_energy = static_cast<floatp_t>(0.0L);
                floatp_t local_energy = static_cast<floatp_t>(0.0L);

                for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
                {
                    const auto position = this->positions[entry_index] - barycenter;

                    const auto rotation_velocity = rotation_axis.cross(position);
                    const auto local_velocity = this->velocities[entry_index] - translation - rotation_velocity;

                    rotation_energy += this->volumes[entry_index] * rotation_velocity.squaredNorm();
                    local_energy += this->volumes[entry_index] * local_velocity.squaredNorm();
                }

                return std::make_pair(rotation_energy, local_energy);
            }
        }
    }
}
