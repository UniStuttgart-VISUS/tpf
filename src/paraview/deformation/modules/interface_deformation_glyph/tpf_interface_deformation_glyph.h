#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

class VTK_EXPORT tpf_interface_deformation_glyph : public vtkMultiBlockDataSetAlgorithm
{
    public:
        static tpf_interface_deformation_glyph* New();
        vtkTypeMacro(tpf_interface_deformation_glyph, vtkMultiBlockDataSetAlgorithm);

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

        vtkSetMacro(BendingDiscResolution, int);
        vtkGetMacro(BendingDiscResolution, int);

        vtkSetMacro(PolygonalResolution, int);
        vtkGetMacro(PolygonalResolution, int);

        vtkSetMacro(BendingStripSize, float);
        vtkGetMacro(BendingStripSize, float);

        vtkSetMacro(BendingSizeScalar, float);
        vtkGetMacro(BendingSizeScalar, float);

        vtkSetMacro(BendingScalar, float);
        vtkGetMacro(BendingScalar, float);

    protected:
        tpf_interface_deformation_glyph();
        ~tpf_interface_deformation_glyph();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

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

        /// Properties of the bending glyph
        int BendingDiscResolution, PolygonalResolution;
        float BendingStripSize, BendingSizeScalar, BendingScalar;
};
