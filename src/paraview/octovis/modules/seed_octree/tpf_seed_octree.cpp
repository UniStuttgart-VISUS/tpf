#include "tpf_seed_octree.h"

#include "tpf/data/tpf_octree.h"
#include "tpf/data/tpf_position.h"

#include "tpf/geometry/tpf_cuboid.h"
#include "tpf/geometry/tpf_point.h"

#include "tpf/log/tpf_log.h"

#include "tpf/mpi/tpf_mpi.h"

#include "tpf/vtk/tpf_data.h"
#include "tpf/vtk/tpf_polydata.h"

#include "vtkDoubleArray.h"
#define FLOAT_TYPE_ARRAY vtkDoubleArray
#ifdef __tpf_no_longlong
#include "vtkLongArray.h"
#define ID_TYPE_ARRAY vtkLongArray
#else
#include "vtkLongLongArray.h"
#define ID_TYPE_ARRAY vtkLongLongArray
#endif

#include "vtkObjectFactory.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

template <typename float_t>
tpf::data::octree<float_t, float_t> load_octree(vtkPointSet* in_octree,
    const std::function<std::string(int)> get_array_name, const bool has_density)
{
    // Get input data
    vtkIdType num_points;

    ID_TYPE_ARRAY* in_paths = nullptr;
    FLOAT_TYPE_ARRAY* in_density = nullptr;

    double x_min, x_max, y_min, y_max, z_min, z_max;

    if (tpf::mpi::get_instance().get_rank() == 0)
    {
        num_points = in_octree->GetPoints()->GetNumberOfPoints();

        in_paths = ID_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(0).c_str()));

        if (has_density)
        {
            in_density = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetPointData()->GetArray(get_array_name(1).c_str()));
        }

        x_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(2).c_str()))->GetValue(0);
        x_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(3).c_str()))->GetValue(0);
        y_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(4).c_str()))->GetValue(0);
        y_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(5).c_str()))->GetValue(0);
        z_min = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(6).c_str()))->GetValue(0);
        z_max = FLOAT_TYPE_ARRAY::SafeDownCast(in_octree->GetFieldData()->GetArray(get_array_name(7).c_str()))->GetValue(0);
    }

#ifdef __tpf_use_mpi
    if (tpf::mpi::get_instance().check_mpi_status())
    {
        tpf::mpi::get_instance().broadcast(num_points, 0);

        if (tpf::mpi::get_instance().get_rank() != 0)
        {
            in_paths = ID_TYPE_ARRAY::New();
            in_paths->SetNumberOfComponents(1);
            in_paths->SetNumberOfTuples(num_points);

            if (has_density)
            {
                in_density = FLOAT_TYPE_ARRAY::New();
                in_density->SetNumberOfComponents(1);
                in_density->SetNumberOfTuples(num_points);
            }
        }

        MPI_Bcast(in_paths->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename ID_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());

        if (has_density)
        {
            MPI_Bcast(in_density->GetVoidPointer(0), num_points, tpf::mpi::mpi_t<typename FLOAT_TYPE_ARRAY::ValueType>::value, 0, tpf::mpi::get_instance().get_comm());
        }

        tpf::mpi::get_instance().broadcast(x_min, 0);
        tpf::mpi::get_instance().broadcast(x_max, 0);
        tpf::mpi::get_instance().broadcast(y_min, 0);
        tpf::mpi::get_instance().broadcast(y_max, 0);
        tpf::mpi::get_instance().broadcast(z_min, 0);
        tpf::mpi::get_instance().broadcast(z_max, 0);
    }
#endif

    // Create octree with root using the bounding box from the data
    const tpf::geometry::point<float_t> min_point(x_min, y_min, z_min);
    const tpf::geometry::point<float_t> max_point(x_max, y_max, z_max);

    auto bounding_box = std::make_shared<tpf::geometry::cuboid<float_t>>(min_point, max_point);
    auto root = std::make_shared<typename tpf::data::octree<float_t, float_t>::node>(std::make_pair(bounding_box, 0.0f));

    tpf::data::octree<float_t, float_t> octree;
    octree.set_root(root);

    // Add nodes
    std::vector<std::pair<std::vector<tpf::data::position_t>, float_t>> node_information;
    node_information.reserve(num_points);

    for (vtkIdType p = 0; p < num_points; ++p)
    {
        // Decode path
        auto path_code = in_paths->GetValue(p);

        std::vector<tpf::data::position_t> path(static_cast<std::size_t>(std::floor(std::log2(path_code) / std::log2(8))), tpf::data::position_t::INVALID);
        std::size_t index = 0;

        while (path_code != 1)
        {
            const auto turn = path_code % 8;

            path[index++] = static_cast<tpf::data::position_t>(((turn & 1) ? 4 : 0) + ((turn & 2) ? 2 : 0) + ((turn & 4) ? 1 : 0));

            path_code /= 8;
        }

        // Add node information
        node_information.push_back(std::make_pair(path, static_cast<float_t>(has_density ? in_density->GetValue(p) : 0.0)));
    }

    // Delete temporary objects
#ifdef __tpf_use_mpi
    if (tpf::mpi::get_instance().get_rank() != 0)
    {
        in_paths->Delete();

        if (has_density)
        {
            in_density->Delete();
        }
    }
#endif

    octree.insert_nodes(node_information.begin(), node_information.end());

    return octree;
}

/// Create seed points at all the leaves
template <typename float_t>
tpf::data::polydata<float_t> seed_at_leaves(const tpf::data::octree<float_t, float_t>& octree)
{
    const auto nodes = octree.get_leaf_nodes();

    tpf::data::polydata<float_t> seed;

    for (const auto& node : nodes)
    {
        seed.insert(std::make_shared<tpf::geometry::point<float_t>>(node->get_value().first->get_center_point()));
    }

    return seed;
}

enum class predicate_t
{
    both, larger, smaller
};

/// Create seed points around the isosurface, defined by the isovalue parameter
template <typename float_t>
tpf::data::polydata<float_t> seed_at_isosurface(const tpf::data::octree<float_t, float_t>& octree,
    const float_t isovalue, const predicate_t predicate = predicate_t::both)
{
    // Find cells, where own and neighboring isovalues are on opposite sides of the input isovalue
    const auto leaves = octree.get_leaf_nodes();

    tpf::data::polydata<float_t> seed;

    for (const auto& leaf : leaves)
    {
        const auto center_isovalue = leaf->get_value().second - isovalue;
        const auto center_position = leaf->get_value().first->get_center_point();

        // Check neighbors' side of the isovalue
        const auto neighbors = octree.get_neighbors(std::static_pointer_cast<typename tpf::data::octree<float_t, float_t>::node>(leaf));

        bool has_positive, has_negative;
        has_positive = has_negative = false;

        for (const auto& neighbor : neighbors)
        {
            has_positive |= (neighbor->get_value().second - isovalue) > 0.0;
            has_negative |= (neighbor->get_value().second - isovalue) < 0.0;
        }

        // If isovalue crosses the threshold, add the position to the seed
        if (((predicate == predicate_t::smaller || predicate == predicate_t::both) && has_positive && center_isovalue < 0.0) ||
            ((predicate == predicate_t::larger || predicate == predicate_t::both) && has_negative && center_isovalue > 0.0))
        {
            seed.insert(std::make_shared<tpf::geometry::point<float_t>>(center_position));
        }
    }

    return seed;
}

vtkStandardNewMacro(tpf_seed_octree);

tpf_seed_octree::tpf_seed_octree()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

tpf_seed_octree::~tpf_seed_octree() {}

int tpf_seed_octree::FillInputPortInformation(int port, vtkInformation* info)
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

    return 0;
}

int tpf_seed_octree::RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *)
{
    return 1;
}

int tpf_seed_octree::RequestData(vtkInformation *request, vtkInformationVector **input_vector, vtkInformationVector *output_vector)
{
    try
    {
        // Get input data and information
        auto in_octree = vtkPointSet::GetData(input_vector[0]);

        const auto get_array_name = [this, &in_octree](const int index) -> std::string
        {
            const auto data = GetInputArrayToProcess(index, in_octree);

            if (data != nullptr)
            {
                return data->GetName();
            }

            throw std::runtime_error(__tpf_error_message("Tried to access non-existing array (ID: ", index, ")."));
        };

        // Sanity check for input
        const auto has_density = GetInputArrayToProcess(1, in_octree) != nullptr;

        if (this->SeedMethod == 1 && !has_density)
        {
            throw std::runtime_error(__tpf_error_message("Iso surface seeding requires density array to be set."));
        }

        // Create octree
        const auto octree = load_octree<float_t>(in_octree, get_array_name, this->SeedMethod == 1 && has_density);

        // Create seed
        tpf::data::polydata<float_t> seed;

        switch (this->SeedMethod)
        {
        case 0:
            seed = seed_at_leaves(octree);
            break;
        case 1:
            seed = seed_at_isosurface(octree, static_cast<float_t>(this->Isovalue), static_cast<predicate_t>(this->Predicate));
            break;
        }

        // Set seed offset and size depending on MPI rank
        std::size_t seed_offset = 0;
        std::size_t seed_size = seed.get_num_objects();

        if (this->SeedOffset >= 0 && this->SeedSize >= 0)
        {
            seed_offset = static_cast<std::size_t>(this->SeedOffset);
            seed_size = static_cast<std::size_t>(this->SeedSize);
        }

#ifdef __tpf_use_mpi
        if (tpf::mpi::get_instance().check_mpi_status())
        {
            const auto num_processes = tpf::mpi::get_instance().get_num_processes();
            const auto rank = tpf::mpi::get_instance().get_rank();

            seed_offset = rank * (seed_size / num_processes);
            seed_size = std::min((rank + 1) * (seed_size / num_processes), seed_size) - seed_offset;
        }
#endif

        tpf::log::info_message(__tpf_info_message("Seed size: ", seed.get_num_objects()), true);
        tpf::log::info_message(__tpf_info_message("Seed range: [", seed_offset, ", ", (seed_offset + seed_size), ")"));

        // Remove seed points to match seed offset and size
        tpf::data::polydata<float_t> reduced_seed;

        for (std::size_t i = seed_offset; i < seed_offset + seed_size; ++i)
        {
            reduced_seed.insert(seed.get_object(i));
        }

        std::swap(seed, reduced_seed);

        // Set output
        auto output = vtkPolyData::GetData(output_vector);

        auto mesh = tpf::vtk::create_mesh(seed.get_geometry());

        output->SetPoints(mesh.points);
        output->SetVerts(mesh.vertices);
    }
    catch (const std::exception& ex)
    {
        tpf::log::error_message(__tpf_nested_error_message(ex.what(), "Execution of plugin 'Seed (Octo-Tiger)' failed."));

        return 0;
    }

    return 1;
}
