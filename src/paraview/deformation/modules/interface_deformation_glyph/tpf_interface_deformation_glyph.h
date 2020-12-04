#pragma once

#include "vtkAlgorithm.h"

class VTK_EXPORT tpf_interface_deformation_glyph : public vtkAlgorithm
{
    public:
        static tpf_interface_deformation_glyph* New();
        vtkTypeMacro(tpf_interface_deformation_glyph, vtkAlgorithm);

        vtkSetMacro(VelocityGlyph, int);
        vtkGetMacro(VelocityGlyph, int);

        vtkSetMacro(StretchingGlyph, int);
        vtkGetMacro(StretchingGlyph, int);

        vtkSetMacro(BendingGlyph, int);
        vtkGetMacro(BendingGlyph, int);

        vtkSetMacro(ArrowSize, int);
        vtkGetMacro(ArrowSize, int);

        vtkSetMacro(ArrowScalar, float);
        vtkGetMacro(ArrowScalar, float);

        vtkSetMacro(ArrowFixedScalar, float);
        vtkGetMacro(ArrowFixedScalar, float);

        vtkSetMacro(Timestep, float);
        vtkGetMacro(Timestep, float);

        vtkSetMacro(ArrowResolution, int);
        vtkGetMacro(ArrowResolution, int);

        vtkSetMacro(ShaftTipRatio, float);
        vtkGetMacro(ShaftTipRatio, float);

        vtkSetMacro(ArrowThickness, float);
        vtkGetMacro(ArrowThickness, float);

        vtkSetMacro(StretchingDiscResolution, int);
        vtkGetMacro(StretchingDiscResolution, int);

        vtkSetMacro(StretchingHoleRadius, float);
        vtkGetMacro(StretchingHoleRadius, float);

        vtkSetMacro(StretchingStripSize, float);
        vtkGetMacro(StretchingStripSize, float);

        vtkSetMacro(StretchingSizeScalar, float);
        vtkGetMacro(StretchingSizeScalar, float);

        vtkSetMacro(BendingDiscResolution, int);
        vtkGetMacro(BendingDiscResolution, int);

        vtkSetMacro(PolynomialResolution, int);
        vtkGetMacro(PolynomialResolution, int);

        vtkSetMacro(BendingStripSize, float);
        vtkGetMacro(BendingStripSize, float);

        vtkSetMacro(BendingSizeScalar, float);
        vtkGetMacro(BendingSizeScalar, float);

        vtkSetMacro(BendingScalar, float);
        vtkGetMacro(BendingScalar, float);

    protected:
        tpf_interface_deformation_glyph();
        ~tpf_interface_deformation_glyph();

        virtual int FillInputPortInformation(int, vtkInformation*) override;
        virtual int FillOutputPortInformation(int, vtkInformation*) override;

        virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

        virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    private:
        tpf_interface_deformation_glyph(const tpf_interface_deformation_glyph&);
        void operator=(const tpf_interface_deformation_glyph&);

        /// Floating point type
        using float_t = float;

        /// Selection of output
        int VelocityGlyph, StretchingGlyph, BendingGlyph;

        /// Custom time step
        float Timestep;

        /// Properties of the velocity glyph
        int ArrowResolution, ArrowSize;
        float ArrowScalar, ArrowFixedScalar, ShaftTipRatio, ArrowThickness;

        /// Properties of the stretching glyph
        int StretchingDiscResolution;
        float StretchingHoleRadius, StretchingStripSize, StretchingSizeScalar;

        /// Properties of the bending glyph
        int BendingDiscResolution, PolynomialResolution;
        float BendingStripSize, BendingSizeScalar, BendingScalar;
};
