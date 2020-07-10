cmake_minimum_required(VERSION 3.12)

pv_module(tpf_flow_field_octree ${PROJECT_NAME} "" tpf_flow_field_octree_target)

target_link_libraries(${tpf_flow_field_octree_target} PRIVATE tpf_import)
