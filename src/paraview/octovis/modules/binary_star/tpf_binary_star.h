#pragma once

#include "vtkAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_binary_star : public vtkAlgorithm
{
    public:
        static tpf_binary_star* New();
        vtkTypeMacro(tpf_binary_star, vtkAlgorithm);

        vtkSetMacro(NumIterations, int);
        vtkGetMacro(NumIterations, int);

        vtkSetMacro(DensityCutoff, double);
        vtkGetMacro(DensityCutoff, double);

    protected:
        tpf_binary_star();
        ~tpf_binary_star();

        virtual int FillInputPortInformation(int, vtkInformation*) override;
        virtual int FillOutputPortInformation(int, vtkInformation*) override;

        virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

        virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    private:
        tpf_binary_star(const tpf_binary_star&);
        void operator=(const tpf_binary_star&);

        /// Using data types
        using float_t = float;

        /// Number of iteration steps
        int NumIterations;

        /// Density cutoff for ignoring low-density cells
        double DensityCutoff;
};
