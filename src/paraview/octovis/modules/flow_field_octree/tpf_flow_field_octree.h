#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"

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

        vtkGetMacro(LocalityMethod, int);
        vtkSetMacro(LocalityMethod, int);

        vtkGetMacro(ForceFixedFrequency, int);
        vtkSetMacro(ForceFixedFrequency, int);

        vtkGetMacro(FrequencyOmega, double);
        vtkSetMacro(FrequencyOmega, double);

        vtkDataArraySelection* GetArraySelection();

        enum class locality_method_t
        {
            none,
            orbit,
            velocity,
            rotation,
            rigid_body
        };

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
        using float_t = double;

        /// Integration parameters
        int Method, NumAdvections, ForceFixedTimeStep;
        double StreamTimeStep;

        /// Locality parameters
        int LocalityMethod, ForceFixedFrequency;
        double FrequencyOmega;

        /// Selection of property arrays
        vtkSmartPointer<vtkDataArraySelection> array_selection;
};
