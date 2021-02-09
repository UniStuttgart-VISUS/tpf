#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_surface_tension : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_surface_tension* New();
        vtkTypeMacro(tpf_surface_tension, vtkRectilinearGridAlgorithm);

        vtkSetMacro(Coefficient, double);
        vtkGetMacro(Coefficient, double);

        vtkSetMacro(Density, double);
        vtkGetMacro(Density, double);

        vtkSetMacro(Timestep, double);
        vtkGetMacro(Timestep, double);

    protected:
        tpf_surface_tension();
        ~tpf_surface_tension();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_surface_tension(const tpf_surface_tension&);
        void operator=(const tpf_surface_tension&);

        /// Floating type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;
        
        /// Surface tension coefficient
        double Coefficient;
        
        /// Density
        double Density;

        /// Time step
        double Timestep;
};
