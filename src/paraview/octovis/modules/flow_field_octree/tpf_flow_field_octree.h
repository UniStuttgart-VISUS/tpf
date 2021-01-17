#pragma once

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT tpf_flow_field_octree : public vtkPolyDataAlgorithm
{
    public:
        static tpf_flow_field_octree* New();
        vtkTypeMacro(tpf_flow_field_octree, vtkPolyDataAlgorithm);

        vtkGetMacro(Method, int);
        vtkSetMacro(Method, int);

        vtkGetMacro(NumAdvections, int);
        vtkSetMacro(NumAdvections, int);

        vtkGetMacro(ForceFixedTimeStep, int);
        vtkSetMacro(ForceFixedTimeStep, int);

        vtkGetMacro(StreamTimeStep, double);
        vtkSetMacro(StreamTimeStep, double);

        vtkGetMacro(ForceFixedFrequency, int);
        vtkSetMacro(ForceFixedFrequency, int);

        vtkGetMacro(FrequencyOmega, double);
        vtkSetMacro(FrequencyOmega, double);
 
    protected:
        tpf_flow_field_octree();
        ~tpf_flow_field_octree();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_flow_field_octree(const tpf_flow_field_octree&);
        void operator=(const tpf_flow_field_octree&);

        /// Using data types
        using float_t = float;

        /// GUI parameters
        int Method, NumAdvections, ForceFixedTimeStep;
        double StreamTimeStep;
        int ForceFixedFrequency;
        double FrequencyOmega;
};
