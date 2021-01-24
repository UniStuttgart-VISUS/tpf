cmake_minimum_required(VERSION 3.12)

pv_module(tpf_binary_star ${PROJECT_NAME} "" tpf_binary_star_target)

target_link_libraries(${tpf_binary_star_target} PRIVATE tpf_import)
