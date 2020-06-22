cmake_minimum_required(VERSION 3.12)

pv_module(tpf_plic ${PROJECT_NAME} "" tpf_plic_target)

target_link_libraries(${tpf_plic_target} PRIVATE tpf_import)
