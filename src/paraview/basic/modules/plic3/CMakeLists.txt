cmake_minimum_required(VERSION 3.12)

pv_module(tpf_plic3 ${PROJECT_NAME} "" tpf_plic3_target)

target_link_libraries(${tpf_plic3_target} PRIVATE tpf_import)
