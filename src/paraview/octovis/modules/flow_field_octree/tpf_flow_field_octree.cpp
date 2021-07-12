#include "tpf_flow_field_octree.h"

#include "tpf_modules/flow_field/tpf_module_flow_field.h"

#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/exception/tpf_not_implemented_exception.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/math/tpf_quaternion.h"

#include "tpf/policies/tpf_interpolatable.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_octree.h"
#include "tpf/vtk/tpf_polydata.h"
#include "tpf/vtk/tpf_time.h"

#include "vtkCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#define FLOAT_TYPE_ARRAY vtkDoubleArray
#ifdef __tpf_no_longlong
#include "vtkLongArray.h"
#define ID_TYPE_ARRAY vtkLongArray
#else
#include "vtkLongLongArray.h"
#define ID_TYPE_ARRAY vtkLongLongArray
#endif
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
    /// Call back class for next time step
    /// <template name="float_t">Floating point type</template>
    /// <template name="point_t">Point type</template>
    template <typename float_t>
    class data_handler : public tpf::modules::flow_field_aux::request_frame_call_back<float_t, tpf::geometry::point<float_t>>
    {
    public:
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="current_timestep">Current time step this plugin was executed in</param>
        /// <param name="timesteps">Available time steps</param>
        /// <param name="request">Original request that will be modified for requesting different time frames</param>
        /// <param name="velocity_alg">Algorithm producing velocity data</param>
        /// <param name="stars_alg">Algorithm producing information about the stars</param>
        /// <param name="get_array_name">Function to retrieve the array names</param>
        /// <param name="locality_method">Method for defining the locality of the stars</param>
        /// <param name="array_selection">Selection of property arrays to interpolate and store at particle positions</param>
        /// <param name="forward">Forward or reverse integration, where true indicates forward, and false reverse direction</param>
        /// <param name="fixed_timestep">Optional fixed timestep</param>
        /// <param name="fixed_frequency">Optional fixed grid angular frequency</param>
        data_handler(const std::size_t current_timestep, const std::vector<float_t>& timesteps, vtkInformation* request,
            vtkAlgorithm* velocity_alg, vtkAlgorithm* stars_alg, std::function<std::string(int, int)> get_array_name,
            const tpf_flow_field_octree::locality_method_t locality_method, vtkDataArraySelection* array_selection, const bool forward,
            std::optional<float_t> fixed_timestep = std::nullopt, std::optional<float_t> fixed_frequency = std::nullopt)
        {
            this->time_offset = this->original_time = current_timestep;
            this->timesteps = timesteps;

            this->request = request;
            this->velocity_alg = velocity_alg;
            this->stars_alg = stars_alg;

            this->get_array_name = get_array_name;

            this->locality_method = locality_method;

            this->array_selection = array_selection;

            this->forward = forward;
            this->fixed_timestep = fixed_timestep;
            this->fixed_frequency = fixed_frequency;
        }

        /// <summary>
        /// Returns the data for the next time step if possible
        /// </summary>
        /// <returns>[Time step delta, velocities, global velocity parts, translation,
        ///  angular velocity, validity, fields to interpolate and store at the particle positions]</returns>
        virtual std::tuple<
            float_t,
            tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)>,
            std::function<bool(const tpf::geometry::point<float_t>&)>,
            std::vector<std::tuple<std::string, std::size_t, tpf::policies::interpolatable_base<tpf::geometry::point<float_t>>*>>> operator()() override
        {
            // Get data
            if ((this->time_offset >= 0 && this->time_offset < this->timesteps.size()) || this->original_time == this->time_offset)
            {
                ID_TYPE_ARRAY* in_classification;
                vtkDataArray* in_star_velocities, * in_star_angular_velocity;

                auto in_octree = vtkPointSet::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));
                auto in_stars = this->stars_alg != nullptr ? vtkPointSet::SafeDownCast(this->stars_alg->GetOutputDataObject(0)) : nullptr;

                auto num_points = in_octree->GetPoints()->GetNumberOfPoints();
                auto points = in_octree->GetPoints();

                auto in_paths = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0, 0).c_str()));

                auto in_velocities = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1, 0).c_str()));

                const auto x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(4, 0).c_str()))->GetValue(0);
                const auto x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(5, 0).c_str()))->GetValue(0);
                const auto y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(6, 0).c_str()))->GetValue(0);
                const auto y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(7, 0).c_str()))->GetValue(0);
                const auto z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(8, 0).c_str()))->GetValue(0);
                const auto z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(9, 0).c_str()))->GetValue(0);

                const auto time = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(10, 0).c_str()))->GetValue(0);

                // Set angular velocity according to selected method
                double angular_frequency;

                if (this->locality_method == tpf_flow_field_octree::locality_method_t::none)
                {
                    angular_frequency = 0.0;
                }
                else if (this->fixed_frequency)
                {
                    angular_frequency = *this->fixed_frequency;
                }
                else
                {
                    angular_frequency = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(11, 0).c_str()))->GetValue(0);
                }

                this->star_position[0].setZero();
                this->star_position[1].setZero();

                // Load data describing translation and rotation of the stars
                if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                {
                    if (in_stars == nullptr)
                    {
                        throw std::runtime_error(__tpf_error_message(
                            "Classification and star information necessary for velocity-based and rigid body locality."));
                    }

                    try
                    {
                        in_classification = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(12, 0).c_str()));

                        if (in_classification == nullptr)
                        {
                            throw std::runtime_error(__tpf_error_message("Classification array is null."));
                        }

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            in_star_velocities = in_stars->GetPointData()->GetArray(get_array_name(13, 1).c_str());

                            if (in_star_velocities == nullptr)
                            {
                                throw std::runtime_error(__tpf_error_message("Classification or star velocity array is null."));
                            }

                            this->star_velocity[0] << in_star_velocities->GetComponent(0, 0), in_star_velocities->GetComponent(0, 1), in_star_velocities->GetComponent(0, 2);
                            this->star_velocity[1] << in_star_velocities->GetComponent(1, 0), in_star_velocities->GetComponent(1, 1), in_star_velocities->GetComponent(1, 2);
                        }

                        if (this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                            this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                        {
                            in_star_angular_velocity = in_stars->GetPointData()->GetArray(get_array_name(14, 1).c_str());

                            if (in_star_angular_velocity == nullptr)
                            {
                                throw std::runtime_error(__tpf_error_message("Star rotation array is null."));
                            }

                            this->angular_velocity[0] << in_star_angular_velocity->GetComponent(0, 0),
                                in_star_angular_velocity->GetComponent(0, 1), in_star_angular_velocity->GetComponent(0, 2);
                            this->angular_velocity[1] << in_star_angular_velocity->GetComponent(1, 0),
                                in_star_angular_velocity->GetComponent(1, 1), in_star_angular_velocity->GetComponent(1, 2);

                            // Get center of mass of the stars
                            std::array<double, 3> tmp_star_pos_1, tmp_star_pos_2;
                            in_stars->GetPoint(0, tmp_star_pos_1.data());
                            in_stars->GetPoint(1, tmp_star_pos_2.data());

                            this->star_position[0] << tmp_star_pos_1[0], tmp_star_pos_1[1], tmp_star_pos_1[2];
                            this->star_position[1] << tmp_star_pos_2[0], tmp_star_pos_2[1], tmp_star_pos_2[2];
                        }
                    }
                    catch (const std::exception& e)
                    {
                        throw std::runtime_error(__tpf_nested_error_message(e.what(),
                            "Classification and star information necessary for velocity-based and rigid body locality."));
                    }
                    catch (...)
                    {
                        throw std::runtime_error(__tpf_error_message(
                            "Classification and star information necessary for velocity-based and rigid body locality."));
                    }
                }

                // Get property arrays
                std::vector<vtkDataArray*> property_arrays;

                for (int index = 0; index < in_octree->GetPointData()->GetNumberOfArrays(); ++index)
                {
                    auto in_array = in_octree->GetPointData()->GetArray(index);

                    if (in_array && in_array->GetName())
                    {
                        if (this->array_selection->ArrayIsEnabled(in_array->GetName()))
                        {
                            property_arrays.push_back(in_array);
                        }
                    }
                }

                // Output dataset information
                std::array<double, 2> vx_range{}, vy_range{}, vz_range{};
                in_velocities->GetValueRange(vx_range.data(), 0);
                in_velocities->GetValueRange(vy_range.data(), 1);
                in_velocities->GetValueRange(vz_range.data(), 2);

                tpf::log::info_message(__tpf_info_message("Number of points: ", num_points));
                tpf::log::info_message(__tpf_info_message("Domain: [", x_min, ", ", y_min, ", ", z_min, "] x [", x_max, ", ", y_max, ", ", z_max, "]"));
                tpf::log::info_message(__tpf_info_message("Velocities: [", vx_range[0], ", ", vy_range[0], ", ", vz_range[0], "] x [", vx_range[1], ", ", vy_range[1], ", ", vz_range[1], "]"));

                // Create octrees
                this->octree = tpf::vtk::get_octree<float_t, float_t, 3>(in_paths,
                    [&in_velocities](vtkIdType p)
                    {
                        return Eigen::Matrix<float_t, 3, 1>(in_velocities->GetComponent(p, 0),
                            in_velocities->GetComponent(p, 1), in_velocities->GetComponent(p, 2));
                    },
                    std::array<double, 6>{ x_min, x_max, y_min, y_max, z_min, z_max });

                if (this->locality_method == tpf_flow_field_octree::locality_method_t::velocity ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rotation ||
                    this->locality_method == tpf_flow_field_octree::locality_method_t::rigid_body)
                {
                    this->octree_classification = tpf::vtk::get_octree<float_t, int, 1>(in_paths,
                        [&in_classification](vtkIdType p)
                        {
                            return static_cast<int>(in_classification->GetValue(p));
                        },
                        std::array<double, 6>{ x_min, x_max, y_min, y_max, z_min, z_max });
                }

                // Create octrees for property arrays
                this->property_octrees.resize(property_arrays.size());

                std::vector<std::tuple<std::string, std::size_t, tpf::policies::interpolatable_base<tpf::geometry::point<float_t>>*>>
                    property_octrees(property_arrays.size());

                for (std::size_t i = 0; i < property_arrays.size(); ++i)
                {
                    if (property_arrays[i]->GetNumberOfComponents() == 3)
                    {
                        this->property_octrees[i] = std::make_shared<tpf::data::octree<float_t, Eigen::Matrix<double, 3, 1>>>(
                            tpf::vtk::get_octree<float_t, double, 3>(in_paths,
                                [&property_arrays, i](vtkIdType p)
                                {
                                    return Eigen::Matrix<double, 3, 1>(property_arrays[i]->GetComponent(p, 0),
                                        property_arrays[i]->GetComponent(p, 1), property_arrays[i]->GetComponent(p, 2));
                                },
                                std::array<double, 6>{ x_min, x_max, y_min, y_max, z_min, z_max }));

                        property_octrees[i] = std::make_tuple(std::string(property_arrays[i]->GetName()), 3, this->property_octrees[i].get());
                    }
                    else if (property_arrays[i]->GetNumberOfComponents() == 1)
                    {
                        this->property_octrees[i] = std::make_shared<tpf::data::octree<float_t, double>>(
                            tpf::vtk::get_octree<float_t, double, 1>(in_paths,
                                [&property_arrays, i](vtkIdType p)
                                {
                                    return property_arrays[i]->GetComponent(p, 0);
                                },
                                std::array<double, 6>{ x_min, x_max, y_min, y_max, z_min, z_max }));

                        property_octrees[i] = std::make_tuple(std::string(property_arrays[i]->GetName()), 1, this->property_octrees[i].get());
                    }
                    else
                    {
                        tpf::log::info_message(__tpf_warning_message("Only property arrays with 1 or 3 components are supported."));
                    }
                }

                // Check if next time step is available
                this->forward ? ++this->time_offset : --this->time_offset;

                if (this->time_offset >= 0 && this->time_offset < this->timesteps.size())
                {
                    // Request update
                    this->request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

                    this->velocity_alg->Update(this->request);
                    this->stars_alg->Update(this->request);
                }

                // Calculate time difference
                auto in_next_grid = vtkPointSet::SafeDownCast(this->velocity_alg->GetOutputDataObject(0));

                const auto next_time = static_cast<float_t>(FLOAT_TYPE_ARRAY::SafeDownCast(
                    in_next_grid->GetFieldData()->GetArray(get_array_name(10, 0).c_str()))->GetValue(0));

                float_t timestep_delta;

                if (this->fixed_timestep)
                {
                    timestep_delta = this->forward ? *this->fixed_timestep : -*this->fixed_timestep;
                }
                else
                {
                    timestep_delta = next_time - time;
                }

                if (timestep_delta == 0.0)
                {
                    timestep_delta = this->forward ? static_cast<float_t>(1.0) : static_cast<float_t>(-1.0);
                }

                tpf::log::info_message(__tpf_info_message("Initial time step size: ", timestep_delta));

                // Create functions for locality information
                auto get_star = [this](const tpf::geometry::point<float_t>& position) -> int
                {
                    return this->octree_classification.find_node(position).first->get_value().second;
                };

                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_translation;
                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_angular_velocity;
                std::function<Eigen::Matrix<float_t, 3, 1>(const tpf::geometry::point<float_t>&)> get_barycenter;
                std::function<bool(const tpf::geometry::point<float_t>&)> is_valid;

                get_barycenter = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_position[get_star(position) - 1]; };

                switch (this->locality_method)
                {
                case tpf_flow_field_octree::locality_method_t::velocity:
                    get_translation = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_velocity[get_star(position) - 1]; };
                    get_angular_velocity = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::rotation:
                    get_translation = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    get_angular_velocity = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[get_star(position) - 1]; };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::rigid_body:
                    get_translation = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->star_velocity[get_star(position) - 1]; };
                    get_angular_velocity = [this, get_star](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[get_star(position) - 1]; };
                    is_valid = [get_star](const tpf::geometry::point<float_t>& position) { return get_star(position) > 0; };

                    break;
                case tpf_flow_field_octree::locality_method_t::orbit:
                case tpf_flow_field_octree::locality_method_t::none:
                default:
                    this->angular_velocity[0] = this->angular_velocity[1] = Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 1.0) * angular_frequency;

                    get_translation = [](const tpf::geometry::point<float_t>& position) { return Eigen::Matrix<float_t, 3, 1>(0.0, 0.0, 0.0); };
                    get_angular_velocity = [this](const tpf::geometry::point<float_t>& position) { return this->angular_velocity[0]; };
                    is_valid = [](const tpf::geometry::point<float_t>& position) { return true; };

                    break;
                }

                // Return [Time step delta, velocities, global velocity parts, translation,
                //   angular velocity, initial translation (dummy), initial angular velocity (dummy), validity]
                return std::make_tuple(timestep_delta,
                    static_cast<tpf::policies::interpolatable<Eigen::Matrix<float_t, 3, 1>, tpf::geometry::point<float_t>>*>(&this->octree),
                    get_translation, get_angular_velocity, get_barycenter, get_translation, get_angular_velocity, is_valid, property_octrees);
            }

            throw std::exception();
        }

        /// <summary>
        /// Reset input algorithms to their beginning state
        /// </summary>
        virtual void reset() override
        {
            // Request update to original time step
            this->time_offset = this->original_time;

            request->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->timesteps[this->time_offset]);

            this->velocity_alg->Update(request);
            this->stars_alg->Update(request);
        }

    private:
        /// Original timestep
        std::size_t original_time;

        /// Last requested timestep
        int time_offset;

        /// Available timesteps
        std::vector<float_t> timesteps;

        /// Request
        vtkInformation* request;

        /// Velocity algorithm
        vtkAlgorithm* velocity_alg;

        /// Star identification algorithm
        vtkAlgorithm* stars_alg;

        /// Function to retrieve array names
        std::function<std::string(int, int)> get_array_name;

        /// Locality definition for the stars
        typename tpf_flow_field_octree::locality_method_t locality_method;

        /// Property arrays
        vtkDataArraySelection* array_selection;

        /// Time direction
        bool forward;

        /// Fixed timestep
        std::optional<float_t> fixed_timestep;

        /// Fixed rotational frequency
        std::optional<float_t> fixed_frequency;

        /// Data of current timestep
        tpf::data::octree<float_t, Eigen::Matrix<float_t, 3, 1>> octree;
        tpf::data::octree<float_t, int> octree_classification;

        std::vector<std::shared_ptr<tpf::policies::interpolatable_base<tpf::geometry::point<float_t>>>> property_octrees;

        /// Position, angular velocity and velocity
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> star_position;
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> angular_velocity;
        std::array<Eigen::Matrix<float_t, 3, 1>, 2> star_velocity;
    };
}

vtkStandardNewMacro(tpf_flow_field_octree);

tpf_flow_field_octree::tpf_flow_field_octree()
{
    this->SetNumberOfInputPorts(3);
    this->SetNumberOfOutputPorts(1);

    this->array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);
}

tpf_flow_field_octree::~tpf_flow_field_octree() {}

vtkDataArraySelection* tpf_flow_field_octree::GetArraySelection()
{
    return this->array_selection;
}

int tpf_flow_field_octree::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }
    if (port == 2)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
        return 1;
    }

    return 0;
}

int tpf_flow_field_octree::RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    return 1;
}

int tpf_flow_field_octree::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data and information
        auto in_octree = vtkPointSet::GetData(input_vector[0]);
        auto in_seed = vtkPolyData::GetData(input_vector[1]);
        auto in_stars = vtkPointSet::GetData(input_vector[2]);

        const auto get_array_name = [this, &in_octree, &in_stars](const int index, const int data_index = 0) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, data_index == 0 ? in_octree : in_stars);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")"));
        };

        const auto current_timestep = tpf::vtk::get_timestep<float_t>(input_vector[0]->GetInformationObject(0), in_octree).second;
        const auto timesteps = tpf::vtk::get_timesteps<float_t>(input_vector[0]->GetInformationObject(0));

        data_handler<float_t> call_back_loader(current_timestep, timesteps, request, GetInputAlgorithm(0, 0), GetInputAlgorithm(2, 0),
            get_array_name, static_cast<locality_method_t>(this->LocalityMethod), this->array_selection, this->Direction == 0,
            this->ForceFixedTimeStep ? std::make_optional(static_cast<float_t>(this->StreamTimeStep)) : std::nullopt,
            this->ForceFixedFrequency ? std::make_optional(static_cast<float_t>(this->FrequencyOmega)) : std::nullopt);

        // Get initial seed
        tpf::data::polydata<float_t> seed = tpf::vtk::get_polydata<float_t>(in_seed);

        // Create and run module for computing integration lines
        using flow_t = tpf::modules::flow_field<float_t, tpf::geometry::point<float_t>>;

        flow_t flow_field;
        tpf::data::polydata<float_t> lines;

        flow_field.set_input(seed);
        flow_field.set_output(lines);

        flow_field.set_parameters(static_cast<tpf::modules::flow_field_aux::method_t>(this->Method),
            static_cast<std::size_t>(this->NumAdvections), tpf::modules::flow_field_aux::time_dependency_t::dynamic, true, true);

        flow_field.set_callbacks(&call_back_loader);

        flow_field.run();

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        tpf::vtk::set_polydata(output, lines,
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Advection)"), tpf::data::topology_t::CELL_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Distribution)"), tpf::data::topology_t::CELL_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Advection)"), tpf::data::topology_t::POINT_DATA },
            tpf::data::data_information<std::size_t, 1>{ std::string("ID (Distribution)"), tpf::data::topology_t::POINT_DATA });

        for (const auto& pd : lines.get_point_data())
        {
            if (pd->get_name().substr(0, 2).compare("ID"))
            {
                tpf::vtk::set_data<double>(output, tpf::data::topology_t::POINT_DATA, pd->get_name(),
                    pd->get_data_dynamic(), pd->get_num_components_dynamic());
            }
        }
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Flow Field (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
