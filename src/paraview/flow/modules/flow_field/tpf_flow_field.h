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

        vtkGetVector2Macro(TimeRange, double);
        vtkSetVector2Macro(TimeRange, double);

        vtkGetMacro(FixedTimeStep, double);
        vtkSetMacro(FixedTimeStep, double);

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

        /// Number of advections for streamlines
        int NumAdvections;

        /// Time range for path- and streaklines
        double TimeRange[2];

        /// Time step size
        double FixedTimeStep;

        /// Dynamic vs. static frame of reference
        int TimeDependency;

        /// Options to keep translational and/or rotational velocity parts
        int KeepTranslation;
        int KeepRotation;

        /// Selection of property arrays
        vtkSmartPointer<vtkDataArraySelection> array_selection;
};
