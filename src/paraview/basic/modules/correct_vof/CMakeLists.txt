cmake_minimum_required(VERSION 3.12)

pv_module(tpf_correct_vof ${PROJECT_NAME} "" tpf_correct_vof_target)

target_link_libraries(${tpf_correct_vof_target} PRIVATE tpf_import)
