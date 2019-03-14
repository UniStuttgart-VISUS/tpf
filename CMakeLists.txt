##################################################################
#
# TPF - Two Phase Flow
#
# Header-only library to fascilitate the development of
# post-processing tools for multiphase flow simulation data,
# integrated in visualization tools, such as ParaView and
# MegaMol.
#
# Copyright 2016-2019
# Visualization Research Center, University of Stuttgart, Germany
#
# Author(s): Alexander Straub
#            Moritz Heinemann
#
##################################################################
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# Create project/solution
project("tpf")

file(STRINGS "version" tpf_version LIMIT_COUNT 1)
message(STATUS "TPF is at version ${tpf_version}")

string(FIND "${tpf_version}" "\." delimiter)
math(EXPR delimiter_plus_one "${delimiter} + 1")

string(SUBSTRING "${tpf_version}" 0 ${delimiter} tpf_version_major)
string(SUBSTRING "${tpf_version}" ${delimiter_plus_one} -1 tpf_version_minor)

# Create header-only target for installation
add_library(tpf INTERFACE)

target_include_directories(tpf INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

install(TARGETS tpf EXPORT tpf-export)

# Copy include directory
install(DIRECTORY "" DESTINATION include/tpf
  PATTERN "*.txt" EXCLUDE
  PATTERN "*.in" EXCLUDE
  PATTERN "*.md" EXCLUDE
  PATTERN ".git" EXCLUDE
  PATTERN "version" EXCLUDE)

# Create config file for including tpf in other projects
configure_file(config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/config.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.cmake DESTINATION share/tpf/cmake RENAME tpf-config.cmake)

# Copy markdown readme files into documentation folder
file(GLOB_RECURSE md_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.md)

foreach(md_file ${md_files})
  get_filename_component(md_file_path ${md_file} DIRECTORY)

  if(md_file_path)
    string(REGEX REPLACE "[/\\]" "_" md_file_name "${md_file_path}")

    install(FILES ${md_file} DESTINATION share/tpf/doc RENAME tpf_${md_file_name}.md)
  else()
    install(FILES ${md_file} DESTINATION share/tpf/doc RENAME tpf.md)
  endif()
endforeach()

# Export project
install(EXPORT tpf-export DESTINATION share/tpf/cmake)
export(TARGETS tpf FILE tpf-export.cmake)

# Create target and source groups for Visual Studio
if(MSVC)
  option(tpf_project_inl_files "Add implementation files to project in order to show them in the IDE" OFF)

  file(GLOB_RECURSE files_and_directories LIST_DIRECTORIES TRUE RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

  set(tpf_files)

  foreach(file_or_directory ${files_and_directories})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${file_or_directory})
      string(REPLACE "/" "\\" filter_name "tpf/${file_or_directory}")

      file(GLOB filter_h_files LIST_DIRECTORIES FALSE "${CMAKE_CURRENT_SOURCE_DIR}/${file_or_directory}/*.h")

      if(tpf_project_inl_files)
        file(GLOB_RECURSE filter_inl_files LIST_DIRECTORIES FALSE "${CMAKE_CURRENT_SOURCE_DIR}/${file_or_directory}/*.inl")
      else()
        set(filter_inl_files)
      endif()

      set(filter_files ${filter_h_files} ${filter_inl_files})
      set(tpf_files ${tpf_files} ${filter_files})

      source_group(${filter_name} FILES ${filter_files})
    endif()
  endforeach()

  add_custom_target(tpf_library SOURCES ${tpf_files})
endif()
