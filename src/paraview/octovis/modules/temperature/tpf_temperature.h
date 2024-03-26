#pragma once

#include "vtkPassInputTypeAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT tpf_temperature : public vtkPassInputTypeAlgorithm
{
    public:
        static tpf_temperature* New();
        vtkTypeMacro(tpf_temperature, vtkPassInputTypeAlgorithm);

    protected:
        tpf_temperature();
        ~tpf_temperature();

        int FillInputPortInformation(int, vtkInformation*);

        int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_temperature(const tpf_temperature&);
        void operator=(const tpf_temperature&);

        /// Using data types
        using float_t = double;
};
