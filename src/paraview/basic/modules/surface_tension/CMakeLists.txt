cmake_minimum_required(VERSION 3.12)

pv_module(tpf_surface_tension ${PROJECT_NAME} "" tpf_surface_tension_target)

target_link_libraries(${tpf_surface_tension_target} PRIVATE tpf_import)
