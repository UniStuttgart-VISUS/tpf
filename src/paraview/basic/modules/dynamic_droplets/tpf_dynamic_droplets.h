#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

class VTK_EXPORT tpf_dynamic_droplets : public vtkMultiBlockDataSetAlgorithm
{
    public:
        static tpf_dynamic_droplets* New();
        vtkTypeMacro(tpf_dynamic_droplets, vtkMultiBlockDataSetAlgorithm);

        vtkSetMacro(NumTimeSteps, int);

        vtkSetMacro(Compute, bool);

        vtkSetMacro(StaticFrameOfReference, bool);

        vtkSetMacro(RibbonSize, double);
        vtkSetMacro(FixRotationAxisSize, bool);
        vtkSetMacro(RotationAxisScale, double);
        vtkSetMacro(RotationAxisTranslation, bool);

    protected:
        tpf_dynamic_droplets();
        ~tpf_dynamic_droplets();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_dynamic_droplets(const tpf_dynamic_droplets&);
        void operator=(const tpf_dynamic_droplets&);

        /// Floating point type
        using float_t = float;

        /// Number of time steps
        int NumTimeSteps;

        /// Option for preventing the costly computation
        bool Compute;

        /// Keep the frame of reference static
        bool StaticFrameOfReference;

        /// Scaling factor for the visualization of ribbons
        double RibbonSize;

        /// Scaling factor and option for the visualization of rotation axes
        bool FixRotationAxisSize;
        double RotationAxisScale;
        bool RotationAxisTranslation;
};
