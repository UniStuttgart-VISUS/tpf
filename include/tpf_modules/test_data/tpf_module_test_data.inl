#include "tpf_module_test_data.h"

#include "tpf_geometry_creator.h"
#include "tpf_orb_creator.h"
#include "tpf_ring_creator.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"
#include "tpf/data/tpf_grid_information.h"

#include "tpf/log/tpf_log.h"

#include <functional>
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
        inline test_data<float_t>::test_data() { }

        template <typename float_t>
        inline std::string test_data<float_t>::get_name() const
        {
            return std::string("Test Data");
        }

        template <typename float_t>
        inline void test_data<float_t>::set_algorithm_output(data::grid<float_t, float_t, 3, 1>& fractions,
            std::optional<std::reference_wrapper<data::grid<float_t, float_t, 3, 3>>> velocities)
        {
            this->fractions = &fractions;
            this->velocities = get_or_default<data::grid<float_t, float_t, 3, 3>>(velocities);
        }

        template <typename float_t>
        inline void test_data<float_t>::set_algorithm_parameters(test_data_aux::geometry_t test_geometry, std::size_t num_cells,
            bool spinning, bool moving, bool expanding, bool rotating, std::optional<float_t> spinning_magnitude,
            std::optional<float_t> moving_magnitude, std::optional<float_t> expanding_magnitude, std::optional<float_t> rotating_magnitude)
        {
            this->test_geometry = test_geometry;

            this->num_cells = num_cells;

            this->spinning = spinning;
            this->moving = moving;
            this->expanding = expanding;
            this->rotating = rotating;

            this->spinning_magnitude = get_or_default(spinning_magnitude, static_cast<float_t>(1.0L));
            this->moving_magnitude = get_or_default(moving_magnitude, static_cast<float_t>(1.0L));
            this->expanding_magnitude = get_or_default(expanding_magnitude, static_cast<float_t>(1.0L));
            this->rotating_magnitude = get_or_default(rotating_magnitude, static_cast<float_t>(1.0L));
        }

        template <typename float_t>
        inline void test_data<float_t>::run_algorithm()
        {
            // Depending on the test case, create different grids
            std::unique_ptr<test_data_aux::geometry_creator<float_t>> creator = nullptr;

            if (this->test_geometry == test_data_aux::ORB)
            {
                create_square_cartesian_grids(this->num_cells);
                creator = std::make_unique<test_data_aux::orb_creator<float_t>>(this->fractions, this->velocities, this->num_cells);
            }
            else if (this->test_geometry == test_data_aux::RING)
            {
                create_square_cartesian_grids(this->num_cells);
                creator = std::make_unique<test_data_aux::ring_creator<float_t>>(this->fractions, this->velocities, this->num_cells);
            }

            // Depending on the test case, add different velocities
            if (creator != nullptr)
            {
                if (this->spinning)
                {
                    creator->add_spinning(this->spinning_magnitude);
                }

                if (this->moving)
                {
                    creator->add_movement(this->moving_magnitude);
                }

                if (this->expanding)
                {
                    creator->add_expansion(this->expanding_magnitude);
                }

                if (this->rotating)
                {
                    creator->add_rotation(this->rotating_magnitude);
                }
            }
        }

        template <typename float_t>
        inline void test_data<float_t>::create_square_cartesian_grids(const std::size_t num_cells)
        {
            // Setup
            const float_t cell_size = static_cast<float_t>(1.0L) / static_cast<float_t>(num_cells);

            typename data::grid_information<float_t>::array_type cell_coordinates, node_coordinates, cell_sizes;
            data::extent_t extent;

            // Calculate grid information
            for (std::size_t c = 0; c < 3; ++c)
            {
                extent.push_back(std::make_pair(static_cast<std::size_t>(0), num_cells - 1));

                cell_coordinates.push_back(std::vector<float_t>(num_cells));
                node_coordinates.push_back(std::vector<float_t>(num_cells + 1));
                cell_sizes.push_back(std::vector<float_t>(num_cells));

                for (std::size_t i = 0; i < num_cells; ++i)
                {
                    cell_coordinates[c][i] = (i + static_cast<float_t>(0.5L)) * cell_size;
                    node_coordinates[c][i] = i * cell_size;
                    cell_sizes[c][i] = cell_size;
                }

                node_coordinates[c][num_cells] = num_cells * cell_size;
            }

            // Create and return empty grids
            *this->fractions = data::grid<float_t, float_t, 3, 1>("VOF", extent, cell_coordinates, node_coordinates, cell_sizes);
            *this->velocities = data::grid<float_t, float_t, 3, 3>("Velocities", extent, cell_coordinates, node_coordinates, cell_sizes);
        }
    }
}
