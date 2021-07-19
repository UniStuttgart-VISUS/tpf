#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_correct_vof : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_correct_vof* New();
        vtkTypeMacro(tpf_correct_vof, vtkRectilinearGridAlgorithm);

    protected:
        tpf_correct_vof();
        ~tpf_correct_vof();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_correct_vof(const tpf_correct_vof&);
        void operator=(const tpf_correct_vof&);

        /// Using data types
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
};
