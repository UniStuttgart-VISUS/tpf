#pragma once

#include "vtkBitArray.h"
#include "vtkCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkShortArray.h"
#include "vtkSignedCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"

namespace tpf
{
    namespace vtk
    {
        template <typename std_t> struct vtk_array { };
        template <> struct vtk_array<char> { typedef vtkCharArray type; };
        template <> struct vtk_array<double> { typedef vtkDoubleArray type; };
        template <> struct vtk_array<float> { typedef vtkFloatArray type; };
        template <> struct vtk_array<int> { typedef vtkIntArray type; };
        template <> struct vtk_array<long> { typedef vtkLongArray type; };
        template <> struct vtk_array<long long> { typedef vtkLongLongArray type; };
        template <> struct vtk_array<short> { typedef vtkShortArray type; };
        template <> struct vtk_array<signed char> { typedef vtkSignedCharArray type; };
        template <> struct vtk_array<unsigned char> { typedef vtkUnsignedCharArray type; };
        template <> struct vtk_array<unsigned int> { typedef vtkUnsignedIntArray type; };
        template <> struct vtk_array<unsigned long> { typedef vtkUnsignedLongArray type; };
        template <> struct vtk_array<unsigned long long> { typedef vtkUnsignedLongLongArray type; };
        template <> struct vtk_array<unsigned short> { typedef vtkUnsignedShortArray type; };

        template <typename vtk_t> struct std_type { };
        template <> struct std_type<vtkBitArray> { typedef int type; };
        template <> struct std_type<vtkCharArray> { typedef char type; };
        template <> struct std_type<vtkDoubleArray> { typedef double type; };
        template <> struct std_type<vtkFloatArray> { typedef float type; };
        template <> struct std_type<vtkIdTypeArray> { typedef vtkIdType type; };
        template <> struct std_type<vtkIntArray> { typedef int type; };
        template <> struct std_type<vtkLongArray> { typedef long type; };
        template <> struct std_type<vtkLongLongArray> { typedef long long type; };
        template <> struct std_type<vtkShortArray> { typedef short type; };
        template <> struct std_type<vtkSignedCharArray> { typedef signed char type; };
        template <> struct std_type<vtkUnsignedCharArray> { typedef unsigned char type; };
        template <> struct std_type<vtkUnsignedIntArray> { typedef unsigned int type; };
        template <> struct std_type<vtkUnsignedLongArray> { typedef unsigned long type; };
        template <> struct std_type<vtkUnsignedLongLongArray> { typedef unsigned long long type; };
        template <> struct std_type<vtkUnsignedShortArray> { typedef unsigned short type; };
    }
}