#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_plic3 : public vtkPolyDataAlgorithm
{
public:
    static tpf_plic3* New();
    vtkTypeMacro(tpf_plic3, vtkPolyDataAlgorithm);

    vtkGetMacro(Epsilon, double);
    vtkSetMacro(Epsilon, double);

    vtkGetMacro(ErrorMarginPlic, double);
    vtkSetMacro(ErrorMarginPlic, double);

    vtkGetMacro(ErrorMarginPlic3, double);
    vtkSetMacro(ErrorMarginPlic3, double);

    vtkGetMacro(NumIterationsPlic, int);
    vtkSetMacro(NumIterationsPlic, int);

    vtkGetMacro(NumIterationsPlic3, int);
    vtkSetMacro(NumIterationsPlic3, int);

    vtkGetMacro(Perturbation, double);
    vtkSetMacro(Perturbation, double);

    vtkGetMacro(HideErrorCubes, int);
    vtkSetMacro(HideErrorCubes, int);

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

    /// Epsilon
    double Epsilon;

    /// Error margin plic
    double ErrorMarginPlic;

    /// Error margin plic3
    double ErrorMarginPlic3;

    /// Number of iterations for PLIC calculation
    int NumIterationsPlic;

    /// Number of iterations for PLIC3 calculation
    int NumIterationsPlic3;

    /// Perturbation
    double Perturbation;

    // Hide error
    int HideErrorCubes;
};
