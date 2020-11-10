cmake_minimum_required(VERSION 3.12)

pv_module(tpf_interface_deformation ${PROJECT_NAME} "" tpf_interface_deformation_target)

target_link_libraries(${tpf_interface_deformation_target} PRIVATE tpf_import)
