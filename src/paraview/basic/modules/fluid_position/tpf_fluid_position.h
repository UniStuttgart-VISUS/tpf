#pragma once

#include "vtkDataObjectAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_fluid_position : public vtkDataObjectAlgorithm
{
public:
    static tpf_fluid_position* New();
    vtkTypeMacro(tpf_fluid_position, vtkDataObjectAlgorithm);

    vtkSetMacro(PositionType, int);
    vtkGetMacro(PositionType, int);

protected:
    tpf_fluid_position();
    ~tpf_fluid_position();

    virtual int FillInputPortInformation(int, vtkInformation*) override;
    virtual int FillOutputPortInformation(int, vtkInformation*) override;

    virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
    virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
    tpf_fluid_position(const tpf_fluid_position&);
    void operator=(const tpf_fluid_position&);

    /// Data type
    using float_t = float;

    /// Number of ghost levels
    std::size_t num_ghost_levels;

    /// Type of positions
    int PositionType;
};
