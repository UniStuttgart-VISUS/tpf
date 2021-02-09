#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "tpf/module/tpf_generic_module_base.h"

#include <memory>
#include <vector>

class VTK_EXPORT tpf_fs3d_reader : public vtkRectilinearGridAlgorithm
{
    public:
        static tpf_fs3d_reader* New();
        vtkTypeMacro(tpf_fs3d_reader, vtkRectilinearGridAlgorithm);

        vtkSetStringMacro(FileName);
        vtkGetStringMacro(FileName);

    protected:
        tpf_fs3d_reader();
        ~tpf_fs3d_reader();

        int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    private:
        tpf_fs3d_reader(const tpf_fs3d_reader&);
        void operator=(const tpf_fs3d_reader&);

        /// Using data types
        using float_t = float;

        /// <summary>
        /// Find closest time step
        /// </summary>
        /// <param name="time">Time</param>
        /// <returns>Time step</returns>
        std::size_t find_timestep(double time);

        /// Module
        std::shared_ptr<tpf::modules::generic_module_base> fs3d_reader;

        /// File name
        char* FileName;

        /// Time steps
        std::vector<double> timesteps;

        /// Number of components
        std::size_t num_components;
};
