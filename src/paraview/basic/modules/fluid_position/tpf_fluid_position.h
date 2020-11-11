#pragma once

#include "vtkAlgorithm.h"

class VTK_EXPORT tpf_fluid_position : public vtkAlgorithm
{
public:
    static tpf_fluid_position* New();
    vtkTypeMacro(tpf_fluid_position, vtkAlgorithm);

    vtkSetMacro(PositionType, int);
    vtkGetMacro(PositionType, int);

protected:
    tpf_fluid_position();
    ~tpf_fluid_position();

    virtual int FillInputPortInformation(int, vtkInformation*) override;
    virtual int FillOutputPortInformation(int, vtkInformation*) override;

    virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
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
