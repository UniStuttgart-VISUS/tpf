#pragma once

#include "vtkRectilinearGridAlgorithm.h"

class VTK_EXPORT tpf_interface_stretching : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_interface_stretching* New();
        vtkTypeMacro(tpf_interface_stretching, vtkRectilinearGridAlgorithm);

    protected:
        tpf_interface_stretching();
        ~tpf_interface_stretching();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_interface_stretching(const tpf_interface_stretching&);
        void operator=(const tpf_interface_stretching&);

        /// Floating point type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
};
