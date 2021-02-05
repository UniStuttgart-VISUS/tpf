#include "tpf_module_fs3d_reader.h"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace tpf
{
    namespace modules
    {
        template <typename float_t>
        inline fs3d_reader<float_t>::fs3d_reader() : fractions(nullptr), velocities(nullptr), timestep(0) { }

        template <typename float_t>
        inline std::string fs3d_reader<float_t>::get_name() const
        {
            return std::string("FS3D Reader");
        }

        template <typename float_t>
        inline void fs3d_reader<float_t>::set_algorithm_output(fs3d_reader_aux::scalar_or_vector<float_t> output)
        {
            if (output.tag == fs3d_reader_aux::scalar_or_vector<float_t>::SCALAR)
            {
                this->fractions = output.scalar_field;
            }
            else
            {
                this->velocities = output.vector_field;
            }
        }

        template <typename float_t>
        inline void fs3d_reader<float_t>::set_algorithm_parameters(const std::string& file_name)
        {
            this->file_name = file_name;

            // Extract grid file name and check for its existance
            this->grid_file_name = this->file_name;

            auto last_path_separator = this->grid_file_name.find_last_of("/\\");
            this->grid_file_name = this->grid_file_name.substr(0, last_path_separator + 1);
            this->grid_file_name += "netz_xyz";

            std::ifstream file;
            file.open(this->grid_file_name + ".dat");

            if (file.good())
            {
                this->grid_file_name += ".dat";
            }
            else
            {
                file.open(this->grid_file_name + ".bin");

                if (file.good())
                {
                    this->grid_file_name += ".bin";
                }
                else
                {
                    throw std::runtime_error(__tpf_error_message("Grid file not found or inaccessible."));
                }
            }

            file.close();

            this->reader = std::make_unique<fs3d_reader_aux::fs3d_data_reader<float_t, float_t>>(this->file_name, this->grid_file_name);
        }

        template <typename float_t>
        inline void fs3d_reader<float_t>::set_run_parameters(const std::size_t timestep, const data::extent_t& extent)
        {
            if (this->file_name == "")
            {
                throw std::runtime_error(__tpf_error_message("Parameters must be assigned before assigning run parameters."));
            }

            this->timestep = timestep;
            this->extent = extent;

            this->reader->set_subextent(extent);
        }

        template <typename float_t>
        inline data_information fs3d_reader<float_t>::provide_information()
        {
            if (this->file_name == "")
            {
                throw std::runtime_error(__tpf_error_message("Module provide information started before assigning parameters."));
            }

            // Read FS3D files for extent information
            return this->reader->get_information();
        }

        template <typename float_t>
        inline void fs3d_reader<float_t>::run_algorithm()
        {
            // Read FS3D files
            if (this->fractions)
            {
                *this->fractions = std::move(this->reader->template read_timestep<1>(this->timestep));
            }

            if (this->velocities)
            {
                *this->velocities = std::move(this->reader->template read_timestep<3>(this->timestep));
            }
        }
    }
}
