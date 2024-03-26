cmake_minimum_required(VERSION 3.12)

pv_module(tpf_temperature ${PROJECT_NAME} "" tpf_temperature_target)

include(CheckTypeSize)
check_type_size("long long" longlong BUILTIN_TYPES_ONLY LANGUAGE CXX)
if (NOT HAVE_longlong)
  target_compile_definitions(${tpf_temperature_target} PUBLIC __tpf_no_longlong)
endif()

target_link_libraries(${tpf_temperature_target} PRIVATE tpf_import)
