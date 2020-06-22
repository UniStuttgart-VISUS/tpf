#include "tpf_module_fs3d_writer.h"

#include "tpf/data/tpf_data_information.h"
#include "tpf/data/tpf_grid.h"

#include "tpf/log/tpf_log.h"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace tpf
{
    namespace modules
    {
        template <typename float_t, int Components>
        inline std::size_t fs3d_writer<float_t, Components>::get_num_required_ghost_levels()
        {
            return 0;
        }

        template <typename float_t, int Components>
        inline fs3d_writer<float_t, Components>::fs3d_writer() : file_name(""), run_parameter_set(false) { }

        template <typename float_t, int Components>
        inline std::string fs3d_writer<float_t, Components>::get_name() const
        {
            return std::string("FS3D Writer");
        }

        template <typename float_t, int Components>
        inline void fs3d_writer<float_t, Components>::set_algorithm_input(const data::grid<float_t, float_t, 3, Components>& grid)
        {
            this->grid = &grid;
        }

        template <typename float_t, int Components>
        inline void fs3d_writer<float_t, Components>::set_algorithm_parameters(const std::string& file_name)
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
                this->grid_file_name += ".bin";
            }

            file.close();

            this->writer = std::make_unique<fs3d_writer_aux::fs3d_data_writer<float_t, float_t>>(this->file_name, this->grid_file_name);
        }

        template <typename float_t, int Components>
        inline void fs3d_writer<float_t, Components>::set_run_parameters(const std::size_t timestep, const data::extent_t& global_extent,
            const double time, const std::string name, const std::string unit, const std::string grid_unit)
        {
            if (this->file_name == "")
            {
                throw std::runtime_error(__tpf_error_message("Parameters must be assigned before assigning run parameters."));
            }

            this->timestep = timestep;
            this->global_extent = global_extent;

            this->time = time;
            this->name = name;
            this->unit = unit;
            this->grid_unit = grid_unit;

            this->writer->set_global_extent(global_extent);

            this->run_parameter_set = true;
        }

        template <typename float_t, int Components>
        inline void fs3d_writer<float_t, Components>::run_algorithm()
        {
            if (!this->run_parameter_set)
            {
                throw std::runtime_error(__tpf_error_message("Module run started before assigning run parameters."));
            }

            // Write FS3D files
            this->writer->write_timestep(this->timestep, this->time, this->name, this->unit, this->grid_unit, *grid);
        }
    }
}
