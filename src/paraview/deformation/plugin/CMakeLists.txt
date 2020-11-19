cmake_minimum_required(VERSION 3.12)

set(parameters)

pv_plugin(${PROJECT_NAME}
  "tpf_interface_bending;tpf_interface_deformation;tpf_interface_deformation_glyph;tpf_interface_stretching"
  "${parameters}"
)
