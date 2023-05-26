cmake_minimum_required(VERSION 3.12)

pv_module(tpf_binary_star ${PROJECT_NAME} "" tpf_binary_star_target)

include(CheckTypeSize)
check_type_size("long long" longlong BUILTIN_TYPES_ONLY LANGUAGE CXX)
if (NOT HAVE_longlong)
  target_compile_definitions(${tpf_binary_star_target} PUBLIC __tpf_no_longlong)
endif()

target_link_libraries(${tpf_binary_star_target} PRIVATE tpf_import)
