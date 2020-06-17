cmake_minimum_required(VERSION 3.12)

# Create ParaView plug-ins
option(tpf_paraview_plugins "Create ParaView plug-ins" OFF)

if(tpf_paraview_plugins)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

  # Create imported interface target for TPF
  include(${CMAKE_BINARY_DIR}/config.cmake)

  add_library(tpf_import INTERFACE IMPORTED)
  set_target_properties(tpf_import PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include;${CGAL_INCLUDE_DIRS};${CGAL_3RD_PARTY_INCLUDE_DIRS};${optional_includes}"
    INTERFACE_COMPILE_OPTIONS "${optional_compile_options}"
    INTERFACE_COMPILE_DEFINITIONS "${optional_compile_definitions}"
    INTERFACE_LINK_LIBRARIES "Eigen3::Eigen;CGAL::CGAL;${optional_link_libraries}"
  )

  # Set build type
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()

  # Find ParaView
  include(GNUInstallDirs)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  set(BUILD_SHARED_LIBS ON)

  find_package(ParaView REQUIRED)

  if(ParaView_VERSION VERSION_LESS 5.7)
    include(${PARAVIEW_USE_FILE})
    include(ParaViewPlugins)
  endif()

  # Get project functions
  include(projects.cmake)

  # Add all subdirectories
  file(GLOB children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

  foreach(child ${children})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child})
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
        message(STATUS "Adding plugin directory '${child}'...")
        add_subdirectory(${child})
      endif()
    endif()
  endforeach()
endif()
