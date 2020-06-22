cmake_minimum_required(VERSION 3.12)

pv_module(tpf_fs3d_writer ${PROJECT_NAME} "" tpf_fs3d_writer_target)

target_link_libraries(${tpf_fs3d_writer_target} PRIVATE tpf_import)
