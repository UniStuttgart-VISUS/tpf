#pragma once

#include "tpf/geometry/tpf_point.h"

#include <vector>

namespace tpf
{
    namespace modules
    {
        namespace flow_field_aux
        {
            /// <summary>
            /// Class for storing particle seeds
            /// </summary>
            /// <template name="floatp_t">Floating point type</template>
            template <typename floatp_t>
            class particle_seed
            {
            public:
                /// <summary>
                /// Use standard constructor
                /// </summary>
                particle_seed() = default;

                /// <summary>
                /// Disallow use of copy constructor
                /// </summary>
                particle_seed(const particle_seed&) = delete;

                /// <summary>
                /// Move constructor
                /// </summary>
                /// <param name="move">Source object</param>
                particle_seed(particle_seed&& move);

                /// <summary>
                /// Get the number of seeds
                /// </summary>
                /// <returns>Number of seeds</returns>
                std::size_t get_size() const;

                /// <summary>
                /// Return the specified seed point
                /// </summary>
                /// <param name="index">Index of the seed</param>
                /// <returns>Seed point</returns>
                const tpf::geometry::point<floatp_t>& get_seed(const std::size_t index) const;

                /// <summary>
                /// Return the specified seed point
                /// </summary>
                /// <param name="index">Index of the seed</param>
                /// <returns>Seed point</returns>
                tpf::geometry::point<floatp_t>& get_seed(const std::size_t index);

                /// <summary>
                /// Initialize seed
                /// </summary>
                /// <template name="forward_iterator_t">Forward iterator type</template>
                /// <param name="begin">First item to insert</param>
                /// <param name="end">Past the last item to insert</param>
                template <typename forward_iterator_t>
                void insert(forward_iterator_t begin, forward_iterator_t end);

            protected:
                /// Seeds
                std::vector<tpf::geometry::point<floatp_t>> seed;
            };
        }
    }
}

#include "tpf_particle_seed.inl"
