cmake_minimum_required(VERSION 3.12)

pv_module(tpf_interface_stretching ${PROJECT_NAME} "" tpf_interface_stretching_target)

target_link_libraries(${tpf_interface_stretching_target} PRIVATE tpf_import)
