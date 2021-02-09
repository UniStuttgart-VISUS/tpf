cmake_minimum_required(VERSION 3.12)

# Add all subdirectories
file(GLOB children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
      message(STATUS "  ...module '${child}'")
      add_subdirectory(${child})
    endif()
  endif()
endforeach()
