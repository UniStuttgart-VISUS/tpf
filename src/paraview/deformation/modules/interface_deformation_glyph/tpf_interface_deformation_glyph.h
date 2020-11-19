#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

class VTK_EXPORT tpf_interface_deformation_glyph : public vtkMultiBlockDataSetAlgorithm
{
    public:
        static tpf_interface_deformation_glyph* New();
        vtkTypeMacro(tpf_interface_deformation_glyph, vtkMultiBlockDataSetAlgorithm);

        vtkSetMacro(Timestep, float);
        vtkGetMacro(Timestep, float);

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

        /// Custom time step
        float Timestep;
};
