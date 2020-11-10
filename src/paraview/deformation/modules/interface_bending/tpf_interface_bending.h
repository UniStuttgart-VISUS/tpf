#pragma once

#include "vtkRectilinearGridAlgorithm.h"

class VTK_EXPORT tpf_interface_bending : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_interface_bending* New();
        vtkTypeMacro(tpf_interface_bending, vtkRectilinearGridAlgorithm);
 
    protected:
        tpf_interface_bending();
        ~tpf_interface_bending();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_interface_bending(const tpf_interface_bending&);
        void operator=(const tpf_interface_bending&);

        /// Floating point type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
};
