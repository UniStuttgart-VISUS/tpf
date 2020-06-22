#include "tpf_geometry_creator.h"

#include "tpf/data/tpf_grid.h"

#include <cstddef>

namespace tpf
{
    namespace modules
    {
        namespace test_data_aux
        {
            template <typename floatp_t>
            inline geometry_creator<floatp_t>::geometry_creator(tpf::data::grid<floatp_t, floatp_t, 3, 1>* fractions,
                tpf::data::grid<floatp_t, floatp_t, 3, 3>* velocities, std::size_t num_cells)
                : fractions(fractions), positions(*velocities), velocities(velocities), num_cells(num_cells) { }
        }
    }
}
