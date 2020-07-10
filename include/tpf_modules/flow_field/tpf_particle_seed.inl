#include "tpf_particle_seed.h"

#include "tpf/geometry/tpf_point.h"

#include <utility>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            template <typename floatp_t>
            inline particle_seed<floatp_t>::particle_seed(particle_seed&& move)
            {
                this->seed = std::move(move.seed);
            }

            template <typename floatp_t>
            inline std::size_t particle_seed<floatp_t>::get_size() const
            {
                return this->seed.size();
            }
            
            template <typename floatp_t>
            inline const tpf::geometry::point<floatp_t>& particle_seed<floatp_t>::get_seed(const std::size_t index) const
            {
                return this->seed[index];
            }
            
            template <typename floatp_t>
            inline tpf::geometry::point<floatp_t>& particle_seed<floatp_t>::get_seed(const std::size_t index)
            {
                return this->seed[index];
            }

            template <typename floatp_t>
            template <typename forward_iterator_t>
            inline void particle_seed<floatp_t>::insert(forward_iterator_t begin, forward_iterator_t end)
            {
                this->seed.insert(this->seed.end(), begin, end);
            }
        }
    }
}
