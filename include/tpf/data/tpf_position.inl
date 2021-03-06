#include "tpf_position.h"

#include <utility>

namespace tpf
{
    namespace data
    {
        inline std::pair<position_t, neighbor_t> operator+(const position_t position, const neighbor_t direction)
        {
            if (position == position_t::INVALID)
            {
                return std::make_pair(position_t::INVALID, neighbor_t::NONE);
            }

            neighbor_t problem_case = neighbor_t::NONE;
            position_t new_position = position;

            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::LEFT))
            {
                if (position == position_t::FRONT_BOTTOM_LEFT || position == position_t::FRONT_TOP_LEFT ||
                    position == position_t::BACK_BOTTOM_LEFT || position == position_t::BACK_TOP_LEFT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::LEFT));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) - 1);
                }
            }
            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::RIGHT))
            {
                if (position == position_t::FRONT_BOTTOM_RIGHT || position == position_t::FRONT_TOP_RIGHT ||
                    position == position_t::BACK_BOTTOM_RIGHT || position == position_t::BACK_TOP_RIGHT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::RIGHT));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) + 1);
                }
            }
            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::BOTTOM))
            {
                if (position == position_t::FRONT_BOTTOM_LEFT || position == position_t::FRONT_BOTTOM_RIGHT ||
                    position == position_t::BACK_BOTTOM_LEFT || position == position_t::BACK_BOTTOM_RIGHT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::BOTTOM));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) - 2);
                }
            }
            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::TOP))
            {
                if (position == position_t::FRONT_TOP_LEFT || position == position_t::FRONT_TOP_RIGHT ||
                    position == position_t::BACK_TOP_LEFT || position == position_t::BACK_TOP_RIGHT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::TOP));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) + 2);
                }
            }
            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::FRONT))
            {
                if (position == position_t::FRONT_BOTTOM_LEFT || position == position_t::FRONT_BOTTOM_RIGHT ||
                    position == position_t::FRONT_TOP_LEFT || position == position_t::FRONT_TOP_RIGHT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::FRONT));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) - 4);
                }
            }
            if (static_cast<int>(direction) & static_cast<int>(neighbor_t::BACK))
            {
                if (position == position_t::BACK_BOTTOM_LEFT || position == position_t::BACK_BOTTOM_RIGHT ||
                    position == position_t::BACK_TOP_LEFT || position == position_t::BACK_TOP_RIGHT)
                {
                    problem_case = static_cast<neighbor_t>(static_cast<int>(problem_case) | static_cast<int>(neighbor_t::BACK));
                }
                else
                {
                    new_position = static_cast<position_t>(static_cast<int>(new_position) + 4);
                }
            }

            return std::make_pair(problem_case == neighbor_t::NONE ? new_position : position_t::INVALID, problem_case);
        }

        inline neighbor_t invert(neighbor_t direction)
        {
            // Get inverse by using xor per direction
            if (static_cast<int>(direction) & (static_cast<int>(neighbor_t::LEFT) | static_cast<int>(neighbor_t::RIGHT)))
            {
                direction = static_cast<neighbor_t>(static_cast<int>(direction) ^ (static_cast<int>(neighbor_t::LEFT) | static_cast<int>(neighbor_t::RIGHT)));
            }
            if (static_cast<int>(direction) & (static_cast<int>(neighbor_t::BOTTOM) | static_cast<int>(neighbor_t::TOP)))
            {
                direction = static_cast<neighbor_t>(static_cast<int>(direction) ^ (static_cast<int>(neighbor_t::BOTTOM) | static_cast<int>(neighbor_t::TOP)));
            }
            if (static_cast<int>(direction) & (static_cast<int>(neighbor_t::FRONT) | static_cast<int>(neighbor_t::BACK)))
            {
                direction = static_cast<neighbor_t>(static_cast<int>(direction) ^ (static_cast<int>(neighbor_t::FRONT) | static_cast<int>(neighbor_t::BACK)));
            }

            return direction;
        }
    }
}
