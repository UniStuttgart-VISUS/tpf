#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"

class VTK_EXPORT tpf_flow_field : public vtkPolyDataAlgorithm
{
    public:
        static tpf_flow_field* New();
        vtkTypeMacro(tpf_flow_field, vtkPolyDataAlgorithm);

        vtkGetMacro(Method, int);
        vtkSetMacro(Method, int);

        vtkGetMacro(SeedInCells, int);
        vtkSetMacro(SeedInCells, int);

        vtkGetMacro(SeedCellType, int);
        vtkSetMacro(SeedCellType, int);

        vtkGetMacro(SeedPerCell, int);
        vtkSetMacro(SeedPerCell, int);

        vtkGetMacro(NumAdvections, int);
        vtkSetMacro(NumAdvections, int);

        vtkGetMacro(ForceFixedTimeStep, int);
        vtkSetMacro(ForceFixedTimeStep, int);

        vtkGetMacro(StreamTimeStep, double);
        vtkSetMacro(StreamTimeStep, double);

        vtkGetMacro(TimeStepFromData, int);
        vtkSetMacro(TimeStepFromData, int);

        vtkGetMacro(TimeDependency, int);
        vtkSetMacro(TimeDependency, int);

        vtkGetMacro(KeepTranslation, int);
        vtkSetMacro(KeepTranslation, int);

        vtkGetMacro(KeepRotation, int);
        vtkSetMacro(KeepRotation, int);

        vtkDataArraySelection* GetArraySelection();

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

        /// Seed in cells instead of using the input seed
        int SeedInCells;

        /// Type of cell in which seeds are created
        int SeedCellType;

        /// Number of seed particles per cell
        int SeedPerCell;

        /// Method to use: 0 - streamlines, 1 - streaklines, 2 - pathlines
        int Method;

        /// Number of advections
        int NumAdvections;

        /// Force fixed time step
        int ForceFixedTimeStep;

        /// Time step if fixed
        double StreamTimeStep;

        /// Force fixed time step
        int TimeStepFromData;

        /// Dynamic vs. static frame of reference
        int TimeDependency;

        /// Options to keep translational and/or rotational velocity parts
        int KeepTranslation;
        int KeepRotation;

        /// Selection of property arrays
        vtkSmartPointer<vtkDataArraySelection> array_selection;
};
