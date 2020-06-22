#pragma once

#include "vtkInformation.h"
#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT tpf_plic3 : public vtkPolyDataAlgorithm
{
public:
    static tpf_plic3* New();
    vtkTypeMacro(tpf_plic3, vtkPolyDataAlgorithm);

    vtkGetMacro(PlicIterations, int);
    vtkSetMacro(PlicIterations, int);

    vtkGetMacro(Plic3Iterations, int);
    vtkSetMacro(Plic3Iterations, int);

protected:
    tpf_plic3();
    ~tpf_plic3();

    int FillInputPortInformation(int port, vtkInformation *info) override;

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    tpf_plic3(const tpf_plic3&) = delete;
    void operator=(const tpf_plic3&) = delete;

    /// Using data types
    using float_t = float;

    /// Number of iterations for PLIC calculation
    int PlicIterations;

    /// Number of iterations for PLIC3 calculation
    int Plic3Iterations;
};
