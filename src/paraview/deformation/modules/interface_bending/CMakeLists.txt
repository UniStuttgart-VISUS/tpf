cmake_minimum_required(VERSION 3.12)

pv_module(tpf_interface_bending ${PROJECT_NAME} "" tpf_interface_bending_target)

target_link_libraries(${tpf_interface_bending_target} PRIVATE tpf_import)
