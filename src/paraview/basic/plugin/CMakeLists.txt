cmake_minimum_required(VERSION 3.12)

set(parameters)
if (ParaView_VERSION VERSION_LESS 5.7)
  set(parameters GUI_RESOURCES tpf_plic.qrc tpf_plic3.qrc)
else()
  set(parameters UI_RESOURCES tpf_plic.qrc tpf_plic3.qrc)
endif()

pv_plugin(${PROJECT_NAME}
  "tpf_correct_vof;tpf_droplets;tpf_dynamic_droplets;tpf_fluid_position;tpf_fs3d_reader;tpf_fs3d_writer;tpf_interface_curvature;tpf_interface_gradient;tpf_plic;tpf_plic3;tpf_stl_reader;tpf_surface_tension;tpf_test_data"
  "${parameters}"
)
