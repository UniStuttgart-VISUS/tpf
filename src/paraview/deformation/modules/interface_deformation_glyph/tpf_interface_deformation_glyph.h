#pragma once

#include "vtkDataObjectAlgorithm.h"

class VTK_EXPORT tpf_interface_deformation_glyph : public vtkDataObjectAlgorithm
{
    public:
        static tpf_interface_deformation_glyph* New();
        vtkTypeMacro(tpf_interface_deformation_glyph, vtkDataObjectAlgorithm);

        vtkSetMacro(VelocityGlyph, int);
        vtkGetMacro(VelocityGlyph, int);

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

        vtkSetMacro(StretchingGlyph, int);
        vtkGetMacro(StretchingGlyph, int);

        vtkSetMacro(StretchingExponent, float);
        vtkGetMacro(StretchingExponent, float);

        vtkSetMacro(StretchingSizeScalar, float);
        vtkGetMacro(StretchingSizeScalar, float);

        vtkSetMacro(StretchingHoleRadius, float);
        vtkGetMacro(StretchingHoleRadius, float);

        vtkSetMacro(StretchingOffset, float);
        vtkGetMacro(StretchingOffset, float);

        vtkSetMacro(StretchingDiscResolution, int);
        vtkGetMacro(StretchingDiscResolution, int);

        vtkSetMacro(StretchingDiscBending, float);
        vtkGetMacro(StretchingDiscBending, float);

        vtkSetMacro(StretchingShowStrips, int);
        vtkGetMacro(StretchingShowStrips, int);

        vtkSetMacro(StretchingStripSize, float);
        vtkGetMacro(StretchingStripSize, float);

        vtkSetMacro(StretchingShowReference, int);
        vtkGetMacro(StretchingShowReference, int);

        vtkSetMacro(StretchingReferenceSize, float);
        vtkGetMacro(StretchingReferenceSize, float);

        vtkSetMacro(StretchingZOffset, float);
        vtkGetMacro(StretchingZOffset, float);

        vtkSetMacro(BendingGlyph, int);
        vtkGetMacro(BendingGlyph, int);

        vtkSetMacro(BendingScalar, float);
        vtkGetMacro(BendingScalar, float);

        vtkSetMacro(BendingSizeScalar, float);
        vtkGetMacro(BendingSizeScalar, float);

        vtkSetMacro(BendingOffset, float);
        vtkGetMacro(BendingOffset, float);

        vtkSetMacro(BendingDiscResolution, int);
        vtkGetMacro(BendingDiscResolution, int);

        vtkSetMacro(PolynomialResolution, int);
        vtkGetMacro(PolynomialResolution, int);

        vtkSetMacro(BendingShowStrips, int);
        vtkGetMacro(BendingShowStrips, int);

        vtkSetMacro(BendingStripSize, float);
        vtkGetMacro(BendingStripSize, float);

        vtkSetMacro(BendingZOffset, float);
        vtkGetMacro(BendingZOffset, float);

    protected:
        tpf_interface_deformation_glyph();
        ~tpf_interface_deformation_glyph();

        virtual int FillInputPortInformation(int, vtkInformation*) override;
        virtual int FillOutputPortInformation(int, vtkInformation*) override;

        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    private:
        tpf_interface_deformation_glyph(const tpf_interface_deformation_glyph&);
        void operator=(const tpf_interface_deformation_glyph&);

        /// Floating point type
        using float_t = float;

        /// Properties of the velocity glyph
        int VelocityGlyph;

        int ArrowSize;
        float ArrowScalar;
        float ArrowFixedScalar;
        float Timestep;
        int ArrowResolution;
        float ShaftTipRatio;
        float ArrowThickness;

        /// Properties of the stretching glyph
        int StretchingGlyph;

        float StretchingExponent;
        float StretchingSizeScalar;
        float StretchingHoleRadius;
        float StretchingOffset;
        int StretchingDiscResolution;
        float StretchingDiscBending;

        int StretchingShowStrips;
        float StretchingStripSize;

        int StretchingShowReference;
        float StretchingReferenceSize;
        float StretchingZOffset;

        /// Properties of the bending glyph
        int BendingGlyph;

        float BendingScalar;
        float BendingSizeScalar;
        float BendingOffset;
        int BendingDiscResolution;
        int PolynomialResolution;

        int BendingShowStrips;
        float BendingStripSize;
        float BendingZOffset;
};
