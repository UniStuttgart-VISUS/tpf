#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_fs3d_writer : public vtkRectilinearGridAlgorithm
{
public:
    static tpf_fs3d_writer* New();
    vtkTypeMacro(tpf_fs3d_writer, vtkRectilinearGridAlgorithm);

    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    vtkSetStringMacro(Name);
    vtkGetStringMacro(Name);

    vtkSetStringMacro(Unit);
    vtkGetStringMacro(Unit);

    vtkSetStringMacro(GridUnit);
    vtkGetStringMacro(GridUnit);

protected:
    tpf_fs3d_writer();
    ~tpf_fs3d_writer();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    tpf_fs3d_writer(const tpf_fs3d_writer&);
    void operator=(const tpf_fs3d_writer&);

    /// Using data types
    using float_t = float;

    /// File name
    char* FileName;

    /// Name and units
    char* Name;
    char* Unit;
    char* GridUnit;

    /// Number of ghost levels
    std::size_t num_ghost_levels;
};
