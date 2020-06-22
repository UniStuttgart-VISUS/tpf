#pragma once

#include "vtkRectilinearGridAlgorithm.h"
 
class VTK_EXPORT tpf_test_data : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_test_data* New();
        vtkTypeMacro(tpf_test_data, vtkRectilinearGridAlgorithm);

        vtkSetMacro(Geometry, int);

        vtkSetMacro(NumCells, int);

        vtkSetMacro(Spinning, int);
        vtkSetMacro(Moving, int);
        vtkSetMacro(Expanding, int);
        vtkSetMacro(Rotating, int);

        vtkSetMacro(SpinningMagnitude, double);
        vtkSetMacro(MovingMagnitude, double);
        vtkSetMacro(ExpandingMagnitude, double);
        vtkSetMacro(RotatingMagnitude, double);
 
    protected:
        tpf_test_data();
        ~tpf_test_data();

        int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_test_data(const tpf_test_data&);
        void operator=(const tpf_test_data&);

        /// Floating point type
        using float_t = float;

        /// Chosen geometry
        int Geometry;

        /// Number of cells
        int NumCells;

        /// Chosen velocities
        int Spinning;
        int Moving;
        int Expanding;
        int Rotating;

        /// Velocity magnitudes
        double SpinningMagnitude;
        double MovingMagnitude;
        double ExpandingMagnitude;
        double RotatingMagnitude;
};
