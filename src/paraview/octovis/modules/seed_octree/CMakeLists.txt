cmake_minimum_required(VERSION 3.12)

pv_module(tpf_seed_octree ${PROJECT_NAME} "" tpf_seed_octree_target)

target_link_libraries(${tpf_seed_octree_target} PRIVATE tpf_import)
