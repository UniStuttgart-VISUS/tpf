#pragma once

#include <utility>

namespace tpf
{
    namespace data
    {
        /// <summary>
        /// Position in a 2x2, or in a 2x2x2 stencil
        /// </summary>
        enum class position_t
        {
            // Invalid position as result of neighbor detection
            INVALID = -1,

            // Positions for 2D, for example in quadtrees
            BOTTOM_LEFT = 0,
            BOTTOM_RIGHT = 1,
            TOP_LEFT = 2,
            TOP_RIGHT = 3,

            // Positions for 3D, for example in octrees
            FRONT_BOTTOM_LEFT = 0,
            FRONT_BOTTOM_RIGHT = 1,
            FRONT_TOP_LEFT = 2,
            FRONT_TOP_RIGHT = 3,
            BACK_BOTTOM_LEFT = 4,
            BACK_BOTTOM_RIGHT = 5,
            BACK_TOP_LEFT = 6,
            BACK_TOP_RIGHT = 7
        };

        /// <summary>
        /// Neighbors in a 2x2, or in a 2x2x2 stencil
        /// </summary>
        enum class neighbor_t
        {
            // No step
            NONE = 0,

            // Single steps
            LEFT = 1,
            RIGHT = 2,
            BOTTOM = 4,
            TOP = 8,
            FRONT = 16,
            BACK = 32,

            // Diagonal steps on a 2D plane
            BOTTOM_LEFT = BOTTOM + LEFT,
            BOTTOM_RIGHT = BOTTOM + RIGHT,
            TOP_LEFT = TOP + LEFT,
            TOP_RIGHT = TOP + RIGHT,

            FRONT_LEFT = FRONT + LEFT,
            FRONT_RIGHT = FRONT + RIGHT,
            BACK_LEFT = BACK + LEFT,
            BACK_RIGHT = BACK + RIGHT,

            FRONT_BOTTOM = FRONT + BOTTOM,
            FRONT_TOP = FRONT + TOP,
            BACK_BOTTOM = BACK + BOTTOM,
            BACK_TOP = BACK + TOP,

            // Diagonal steps in a 3D cube
            FRONT_BOTTOM_LEFT = FRONT + BOTTOM + LEFT,
            FRONT_BOTTOM_RIGHT = FRONT + BOTTOM + RIGHT,
            FRONT_TOP_LEFT = FRONT + TOP + LEFT,
            FRONT_TOP_RIGHT = FRONT + TOP + RIGHT,
            BACK_BOTTOM_LEFT = BACK + BOTTOM + LEFT,
            BACK_BOTTOM_RIGHT = BACK + BOTTOM + RIGHT,
            BACK_TOP_LEFT = BACK + TOP + LEFT,
            BACK_TOP_RIGHT = BACK + TOP + RIGHT
        };

        /// <summary>
        /// Operator for movement
        /// </summary>
        /// <param name="position">Position</param>
        /// <param name="direction">Direction</param>
        /// <returns>New position</returns>
        std::pair<position_t, neighbor_t> operator+(position_t position, neighbor_t direction);
        neighbor_t invert(neighbor_t direction);
    }
}

#include "tpf_position.inl"
