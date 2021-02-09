#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_interface_gradient : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_interface_gradient* New();
        vtkTypeMacro(tpf_interface_gradient, vtkRectilinearGridAlgorithm);

    protected:
        tpf_interface_gradient();
        ~tpf_interface_gradient();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_interface_gradient(const tpf_interface_gradient&);
        void operator=(const tpf_interface_gradient&);

        /// Using data types
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
};
