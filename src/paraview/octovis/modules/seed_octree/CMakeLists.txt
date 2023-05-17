cmake_minimum_required(VERSION 3.12)

pv_module(tpf_seed_octree ${PROJECT_NAME} "" tpf_seed_octree_target)

include(CheckTypeSize)
check_type_size("long long" longlong BUILTIN_TYPES_ONLY LANGUAGE CXX)
if (NOT HAVE_longlong)
  target_compile_definitions(${tpf_seed_octree_target} PUBLIC __tpf_no_longlong)
endif()

target_link_libraries(${tpf_seed_octree_target} PRIVATE tpf_import)
