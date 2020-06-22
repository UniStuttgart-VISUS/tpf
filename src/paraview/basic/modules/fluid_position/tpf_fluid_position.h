#pragma once

#include "vtkMultiBlockDataSetAlgorithm.h"

class VTK_EXPORT tpf_fluid_position : public vtkMultiBlockDataSetAlgorithm
{
public:
    static tpf_fluid_position* New();
    vtkTypeMacro(tpf_fluid_position, vtkMultiBlockDataSetAlgorithm);

    vtkSetMacro(PositionType, int);
    vtkGetMacro(PositionType, int);

protected:
    tpf_fluid_position();
    ~tpf_fluid_position();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

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
