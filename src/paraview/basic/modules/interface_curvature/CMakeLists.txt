cmake_minimum_required(VERSION 3.12)

pv_module(tpf_interface_curvature ${PROJECT_NAME} "" tpf_interface_curvature_target)

target_link_libraries(${tpf_interface_curvature_target} PRIVATE tpf_import)
