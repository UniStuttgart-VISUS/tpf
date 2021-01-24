#pragma once

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT tpf_binary_star : public vtkPolyDataAlgorithm
{
    public:
        static tpf_binary_star* New();
        vtkTypeMacro(tpf_binary_star, vtkPolyDataAlgorithm);

        vtkGetMacro(SeedMethod, int);
        vtkSetMacro(SeedMethod, int);

        vtkGetMacro(SeedOffset, int);
        vtkSetMacro(SeedOffset, int);

        vtkGetMacro(SeedSize, int);
        vtkSetMacro(SeedSize, int);

        vtkGetMacro(Isovalue, double);
        vtkSetMacro(Isovalue, double);
 
    protected:
        tpf_binary_star();
        ~tpf_binary_star();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_binary_star(const tpf_binary_star&);
        void operator=(const tpf_binary_star&);

        /// Using data types
        using float_t = float;

        /// GUI parameters
        int SeedMethod, SeedOffset, SeedSize;
        double Isovalue;
};
