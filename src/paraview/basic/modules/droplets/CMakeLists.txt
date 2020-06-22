cmake_minimum_required(VERSION 3.12)

pv_module(tpf_droplets ${PROJECT_NAME} "" tpf_droplets_target)

target_link_libraries(${tpf_droplets_target} PRIVATE tpf_import)
