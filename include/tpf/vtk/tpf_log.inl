#include "tpf_log.h"

#include "../log/tpf_log.h"
#include "../log/tpf_log_base.h"

#ifdef __tpf_use_paraview_output_window
#include "vtkObject.h"
#include "vtkOutputWindow.h"
#endif

#include <functional>
#include <memory>
#include <string>

namespace tpf
{
    namespace vtk
    {
        inline vtk_log::vtk_log(std::function<void(const char*)> func, bool active)
        {
            this->func = func;
            this->active = active;
        }

        inline void operator<<(vtk_log& stream, const std::string& input)
        {
            if (stream.active)
            {
                stream.func(input.c_str());
            }
        }

        inline void vtk_log::flush() const {}

        inline std::shared_ptr<log::log_base> create_log_instance()
        {
#ifdef __tpf_use_paraview_output_window
            if (vtkObject::GetGlobalWarningDisplay())
            {
                static vtk_log err_stream(static_cast<std::function<void(const char*)>>(vtkOutputWindowDisplayErrorText));
                static vtk_log warn_stream(static_cast<std::function<void(const char*)>>(vtkOutputWindowDisplayWarningText));
                static vtk_log info_stream(static_cast<std::function<void(const char*)>>(vtkOutputWindowDisplayDebugText));

                return tpf::log::create_instance(err_stream, warn_stream, info_stream);
            }
            else
            {
#endif
                return tpf::log::create_instance();
#ifdef __tpf_use_paraview_output_window
            }
#endif
        }
    }
}