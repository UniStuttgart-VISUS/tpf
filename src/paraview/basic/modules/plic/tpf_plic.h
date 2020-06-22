#pragma once

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT tpf_plic : public vtkPolyDataAlgorithm
{
public:
    static tpf_plic* New();
    vtkTypeMacro(tpf_plic, vtkPolyDataAlgorithm);

    vtkGetMacro(NumIterations, int);
    vtkSetMacro(NumIterations, int);

    vtkGetMacro(Perturbation, double);
    vtkSetMacro(Perturbation, double);

protected:
    tpf_plic();
    ~tpf_plic();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    tpf_plic(const tpf_plic&);
    void operator=(const tpf_plic&);

    /// Using data types
    using float_t = float;

    /// Number of ghost levels
    std::size_t num_ghost_levels;

    /// Number of iterations
    int NumIterations;

    /// Perturbation
    double Perturbation;
};
