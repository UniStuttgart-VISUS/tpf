cmake_minimum_required(VERSION 3.12)

pv_module(tpf_test_data ${PROJECT_NAME} "" tpf_test_data_target)

target_link_libraries(${tpf_test_data_target} PRIVATE tpf_import)
