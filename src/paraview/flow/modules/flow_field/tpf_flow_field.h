#pragma once

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT tpf_flow_field : public vtkPolyDataAlgorithm
{
    public:
        static tpf_flow_field* New();
        vtkTypeMacro(tpf_flow_field, vtkPolyDataAlgorithm);

        vtkGetMacro(Method, int);
        vtkSetMacro(Method, int);

        vtkGetMacro(SeedPerCell, int);
        vtkSetMacro(SeedPerCell, int);

        vtkGetMacro(NumAdvections, int);
        vtkSetMacro(NumAdvections, int);

        vtkGetMacro(ForceFixedTimeStep, int);
        vtkSetMacro(ForceFixedTimeStep, int);

        vtkGetMacro(StreamTimeStep, double);
        vtkSetMacro(StreamTimeStep, double);

        enum class locality_method_t
        {
            none,
            velocity,
            rotation,
            rigid_body
        };

    protected:
        tpf_flow_field();
        ~tpf_flow_field();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_flow_field(const tpf_flow_field&);
        void operator=(const tpf_flow_field&);

        /// Floating point type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;

        /// Method to use: 0 - streamlines, 1 - streaklines, 2 - pathlines
        int Method;

        /// Number of seed particles per cell
        int SeedPerCell;

        /// Number of advections
        int NumAdvections;

        /// Force fixed time step
        int ForceFixedTimeStep;

        /// Time step if fixed
        double StreamTimeStep;
};
