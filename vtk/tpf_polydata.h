#pragma once

#include "../data/tpf_polydata.h"

#include "vtkPolyData.h"

namespace tpf
{
    namespace vtk
    {
        /// <summary>
        /// Get the poly data object
        /// </summary>
        /// <template name="point_t">Point type</template>
        /// <template name="data_info_t">Data information object type (must be data_information)</template>
        /// <param name="data">VTK poly data</param>
        /// <param name="data_info">Data information objects</param>
        /// <returns>Poly data</returns>
        template <typename point_t, typename... data_info_t>
        data::polydata<point_t> get_polydata(vtkPolyData* data, const data_info_t&... data_info);

        /// <summary>
        /// Set the poly data object
        /// </summary>
        /// <template name="point_t">Point type</template>
        /// <template name="data_info_t">Data information object type (must be data_information)</template>
        /// <param name="data">VTK poly data</param>
        /// <param name="polydata">Poly data</param>
        /// <param name="data_info">Data information objects</param>
        template <typename point_t, typename... data_info_t>
        void set_polydata(vtkPolyData* data, const data::polydata<point_t>& polydata, const data_info_t&... data_info);

        /// <summary>
        /// Append to the poly data object
        /// </summary>
        /// <template name="point_t">Point type</template>
        /// <template name="data_info_t">Data information object type (must be data_information)</template>
        /// <param name="data">VTK poly data</param>
        /// <param name="polydata">Poly data</param>
        /// <param name="data_info">Data information objects</param>
        template <typename point_t, typename... data_info_t>
        void append_polydata(vtkPolyData* data, const data::polydata<point_t>& polydata, const data_info_t&... data_info);
    }
}

#include "tpf_polydata.inl"