#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "tpf/module/tpf_generic_module_base.h"

#include <memory>
 
class VTK_EXPORT tpf_stl_reader : public vtkPolyDataAlgorithm
{
    public:
        static tpf_stl_reader* New();
        vtkTypeMacro(tpf_stl_reader, vtkPolyDataAlgorithm);

        vtkSetStringMacro(FileName);
        vtkGetStringMacro(FileName);

        vtkSetMacro(CalculateNormals, int);
        vtkGetMacro(CalculateNormals, int);
 
    protected:
        tpf_stl_reader();
        ~tpf_stl_reader();

        int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
 
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_stl_reader(const tpf_stl_reader&);
        void operator=(const tpf_stl_reader&);

        /// Using data types
        using float_t = float;

        /// Module
        std::shared_ptr<tpf::modules::generic_module_base> stl_reader;

        /// File name
        char* FileName;

        /// Calculate normals instead of reading them
        int CalculateNormals;

        /// Number of triangles
        std::size_t num_triangles;
};
