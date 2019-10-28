#pragma once

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
            BOTTOM_LEFT = 5,
            BOTTOM_RIGHT = 6,
            TOP_LEFT = 9,
            TOP_RIGHT = 10,
            FRONT_LEFT = 17,
            FRONT_RIGHT = 18,
            FRONT_BOTTOM = 20,
            FRONT_TOP = 24,
            BACK_LEFT = 33,
            BACK_RIGHT = 34,
            BACK_BOTTOM = 36,
            BACK_TOP = 40,

            // Diagonal steps in a 3D cube
            FRONT_BOTTOM_LEFT = 21,
            FRONT_BOTTOM_RIGHT = 22,
            FRONT_TOP_LEFT = 25,
            FRONT_TOP_RIGHT = 26,
            BACK_BOTTOM_LEFT = 37,
            BACK_BOTTOM_RIGHT = 38,
            BACK_TOP_LEFT = 41,
            BACK_TOP_RIGHT = 42
        };

        /// <summary>
        /// Operator for movement
        /// </summary>
        /// <param name="position">Position</param>
        /// <param name="direction">Direction</param>
        /// <returns>New position</returns>
        position_t operator+(position_t position, neighbor_t direction);
        position_t operator-(position_t position, neighbor_t direction);
    }
}

#include "tpf_position.inl"