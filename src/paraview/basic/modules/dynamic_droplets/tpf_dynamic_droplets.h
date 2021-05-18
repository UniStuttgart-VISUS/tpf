#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

class VTK_EXPORT tpf_dynamic_droplets : public vtkMultiBlockDataSetAlgorithm
{
    public:
        static tpf_dynamic_droplets* New();
        vtkTypeMacro(tpf_dynamic_droplets, vtkMultiBlockDataSetAlgorithm);

        vtkSetMacro(NumTimeSteps, int);

        vtkSetMacro(ShowTranslationPaths, bool);
        vtkSetMacro(ShowRotationAxes, bool);
        vtkSetMacro(ShowRotationRibbons, bool);

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

        /// Options for different visualizations
        bool ShowTranslationPaths, ShowRotationAxes, ShowRotationRibbons;

        /// Scaling factor and option for the visualization of rotation axes
        double RotationAxisScale;
        bool RotationAxisTranslation;
};
