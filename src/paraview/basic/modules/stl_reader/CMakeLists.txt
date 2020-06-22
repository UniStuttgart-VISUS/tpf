cmake_minimum_required(VERSION 3.12)

pv_module(tpf_stl_reader ${PROJECT_NAME} "" tpf_stl_reader_target)

target_link_libraries(${tpf_stl_reader_target} PUBLIC tpf_import)
