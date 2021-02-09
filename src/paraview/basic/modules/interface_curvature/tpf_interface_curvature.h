#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_interface_curvature : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_interface_curvature* New();
        vtkTypeMacro(tpf_interface_curvature, vtkRectilinearGridAlgorithm);

    protected:
        tpf_interface_curvature();
        ~tpf_interface_curvature();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_interface_curvature(const tpf_interface_curvature&);
        void operator=(const tpf_interface_curvature&);

        /// Floating point type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
};
