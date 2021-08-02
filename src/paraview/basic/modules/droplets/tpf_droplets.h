#pragma once

#include "vtkDataObjectAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_droplets : public vtkDataObjectAlgorithm
{
    public:
        static tpf_droplets* New();
        vtkTypeMacro(tpf_droplets, vtkDataObjectAlgorithm);

        vtkGetMacro(CalculateTranslation, int);
        vtkSetMacro(CalculateTranslation, int);

        vtkGetMacro(CalculateRotation, int);
        vtkSetMacro(CalculateRotation, int);

        vtkGetMacro(CalculateEnergy, int);
        vtkSetMacro(CalculateEnergy, int);

        vtkGetMacro(CalculateInertia, int);
        vtkSetMacro(CalculateInertia, int);

        vtkGetMacro(Scale, int);
        vtkSetMacro(Scale, int);

        vtkGetMacro(RotationMethod, int);
        vtkSetMacro(RotationMethod, int);

    protected:
        tpf_droplets();
        ~tpf_droplets();

        virtual int FillInputPortInformation(int, vtkInformation*) override;
        virtual int FillOutputPortInformation(int, vtkInformation*) override;

        virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    private:
        tpf_droplets(const tpf_droplets&);
        void operator=(const tpf_droplets&);

        /// Floating point type
        using float_t = float;

        /// Number of ghost levels
        std::size_t num_ghost_levels;

        /// Deciders on what to calculate
        int CalculateTranslation;
        int CalculateRotation;
        int CalculateEnergy;
        int CalculateInertia;

        /// Scaling for inertia
        int Scale;

        /// Method for calculating rotation
        int RotationMethod;
};
