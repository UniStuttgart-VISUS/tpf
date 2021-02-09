#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_seed_octree : public vtkPolyDataAlgorithm
{
    public:
        static tpf_seed_octree* New();
        vtkTypeMacro(tpf_seed_octree, vtkPolyDataAlgorithm);

        vtkGetMacro(SeedMethod, int);
        vtkSetMacro(SeedMethod, int);

        vtkGetMacro(SeedOffset, int);
        vtkSetMacro(SeedOffset, int);

        vtkGetMacro(SeedSize, int);
        vtkSetMacro(SeedSize, int);

        vtkGetMacro(Isovalue, double);
        vtkSetMacro(Isovalue, double);
 
    protected:
        tpf_seed_octree();
        ~tpf_seed_octree();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_seed_octree(const tpf_seed_octree&);
        void operator=(const tpf_seed_octree&);

        /// Using data types
        using float_t = double;

        /// GUI parameters
        int SeedMethod, SeedOffset, SeedSize;
        double Isovalue;
};
