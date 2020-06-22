cmake_minimum_required(VERSION 3.12)

# Create module
pv_module(tpf_fluid_position ${PROJECT_NAME} "" tpf_fluid_position_target)

target_link_libraries(${tpf_fluid_position_target} PRIVATE tpf_import)
