cmake_minimum_required(VERSION 3.12)

pv_module(tpf_interface_gradient ${PROJECT_NAME} "" tpf_interface_gradient_target)

target_link_libraries(${tpf_interface_gradient_target} PRIVATE tpf_import)
