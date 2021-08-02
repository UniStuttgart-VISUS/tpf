#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_binary_star : public vtkPolyDataAlgorithm
{
    public:
        static tpf_binary_star* New();
        vtkTypeMacro(tpf_binary_star, vtkPolyDataAlgorithm);

        vtkSetMacro(NumIterations, int);
        vtkGetMacro(NumIterations, int);

        vtkSetMacro(DensityCutoff, double);
        vtkGetMacro(DensityCutoff, double);

    protected:
        tpf_binary_star();
        ~tpf_binary_star();

        virtual int FillInputPortInformation(int, vtkInformation*) override;

        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    private:
        tpf_binary_star(const tpf_binary_star&);
        void operator=(const tpf_binary_star&);

        /// Using data types
        using float_t = double;

        /// Number of iteration steps
        int NumIterations;

        /// Density cutoff for ignoring low-density cells
        double DensityCutoff;
};
