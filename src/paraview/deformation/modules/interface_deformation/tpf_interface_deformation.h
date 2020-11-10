#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include <algorithm>

class VTK_EXPORT tpf_interface_deformation : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_interface_deformation* New();
        vtkTypeMacro(tpf_interface_deformation, vtkRectilinearGridAlgorithm);

        vtkSetMacro(SurfaceTension, int);
        vtkGetMacro(SurfaceTension, int);

        vtkSetMacro(Coefficient, double);
        vtkGetMacro(Coefficient, double);

        vtkSetMacro(Density, double);
        vtkGetMacro(Density, double);

        vtkSetMacro(Timestep, double);
        vtkGetMacro(Timestep, double);
 
    protected:
        tpf_interface_deformation();
        ~tpf_interface_deformation();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_interface_deformation(const tpf_interface_deformation&);
        void operator=(const tpf_interface_deformation&);

        /// Using data types
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;

        /// Parameters
        int SurfaceTension;
        double Coefficient, Density, Timestep;
};
