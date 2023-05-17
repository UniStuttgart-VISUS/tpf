cmake_minimum_required(VERSION 3.12)

pv_module(tpf_flow_field ${PROJECT_NAME} "" tpf_flow_field_target)

include(CheckTypeSize)
check_type_size("long long" longlong BUILTIN_TYPES_ONLY LANGUAGE CXX)
if (NOT HAVE_longlong)
  target_compile_definitions(${tpf_flow_field_target} PUBLIC __tpf_no_longlong)
endif()

target_link_libraries(${tpf_flow_field_target} PRIVATE tpf_import)
