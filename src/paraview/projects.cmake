cmake_minimum_required(VERSION 3.12)

function(pv_project NAME ENABLED)
  project(${NAME} CXX)

  if (ParaView_VERSION VERSION_LESS 5.7)
    option(PARAVIEW_PLUGIN_ENABLE_${NAME} "Enable the ${NAME} plugin." ${ENABLED})

    if (PARAVIEW_PLUGIN_ENABLE_${NAME})
      # Build modules
      add_subdirectory(modules)

      # Build plugin
      add_subdirectory(plugin)
    endif()
  else()
    # Scan plugin
    paraview_plugin_scan(
      PLUGIN_FILES                "${CMAKE_CURRENT_SOURCE_DIR}/plugin/${PROJECT_NAME}.plugin"
      PROVIDES_PLUGINS            plugins
      ENABLE_BY_DEFAULT           ${ENABLED}
      HIDE_PLUGINS_FROM_CACHE     OFF
    )

    if (PARAVIEW_PLUGIN_ENABLE_${NAME})
      # Scan and build modules
      file(GLOB_RECURSE module_files modules/*.module)

      vtk_module_scan(
        MODULE_FILES                ${module_files}
        PROVIDES_MODULES            modules
        WANT_BY_DEFAULT             ON
        HIDE_MODULES_FROM_CACHE     OFF
      )

      foreach(module ${modules})
        message(STATUS "  ...module '${module}'")
      endforeach()

      vtk_module_build(
        MODULES                     ${modules}
        PACKAGE                     ${PROJECT_NAME}
        INSTALL_HEADERS             OFF
      )

      # Build plugin
      message(STATUS "  ...plugin '${PROJECT_NAME}'")

      paraview_plugin_build(
        PLUGINS                     ${plugins}
      )
    endif()
  endif()
endfunction()

function(pv_plugin NAME MODULES PARAMETERS)
  # Collect modules and corresponding headers
  set(module_names)
  set(module_headers)

  foreach(module ${MODULES})
    list(APPEND module_names ${NAME}::${module})

    get_target_property(_module_headers ${NAME}::${module} VTK_HEADERS)
    list(APPEND module_headers ${_module_headers})
  endforeach()

  # Create plugin
  if (ParaView_VERSION VERSION_LESS 5.7)
    message(STATUS "  ...plugin '${NAME}'")

    add_paraview_plugin(${NAME} 1.0
      SERVER_MANAGER_XML ${NAME}.xml
      SERVER_MANAGER_SOURCES ${module_headers} tpf_main.cpp
      ${PARAMETERS}
    )

    target_link_libraries(${NAME} PRIVATE ${module_names} tpf_import)

    install(TARGETS ${NAME})
    install(FILES ${NAME}.xml DESTINATION share/tpf/xml)
  else()
    paraview_add_plugin(${NAME}
      VERSION                  1.0
      SERVER_MANAGER_XML      ${NAME}.xml
      ${PARAMETERS}
      MODULES                 ${module_names}
      SOURCES                 tpf_main.cpp
    )

    target_link_libraries(${NAME} PRIVATE tpf_import)

    install(FILES ${NAME}.xml DESTINATION share/tpf/xml)
  endif()
endfunction()

function(pv_module NAME PLUGIN SOURCES RESULT_TARGET)
  if (ParaView_VERSION VERSION_LESS 5.7)
    # Create module as static library
    add_library(${PLUGIN}_${NAME} STATIC ${NAME}.h ${NAME}.cpp ${SOURCES})
    add_library(${PLUGIN}::${NAME} ALIAS ${PLUGIN}_${NAME})

    target_link_libraries(${PLUGIN}_${NAME} PRIVATE ${VTK_LIBRARIES})
    target_include_directories(${PLUGIN}_${NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${common_include})

    set_target_properties(${PLUGIN}_${NAME} PROPERTIES VTK_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.h")
    set_target_properties(${PLUGIN}_${NAME} PROPERTIES CXX_STANDARD 17) # Revise warning settings below when updating

    # Silence warnings
    target_compile_definitions(${PLUGIN}_${NAME} PRIVATE
      # Deriving from std::iterator is deprecated in C++17 but used by VTK/ParaView. Revise when updating C++ standard!
      _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
      # std::unary_negate and std::binary_negate are deprecated in C++17 but used by Eigen. Revise when updating C++ standard!
      _SILENCE_CXX17_NEGATORS_DEPRECATION_WARNING
      # Many result_type typedefs are deprecated in C++17 but used by Eigen. Revise when updating C++ standard!
      _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)

    set(${RESULT_TARGET} ${PLUGIN}_${NAME} PARENT_SCOPE)
  else()
    # Create VTK module
    vtk_module_add_module(${PLUGIN}::${NAME} FORCE_STATIC
      SOURCES
        ${NAME}.cpp
        ${SOURCES}
      HEADERS
        ${NAME}.h
    )

    _vtk_module_real_target(_RESULT_TARGET ${PLUGIN}::${NAME})
    set(${RESULT_TARGET} ${_RESULT_TARGET} PARENT_SCOPE)

    target_include_directories(${_RESULT_TARGET} PRIVATE ${common_include})
    vtk_module_set_properties(${PLUGIN}::${NAME} CXX_STANDARD 17) # Revise warning settings below when updating

    # Silence warnings
    vtk_module_definitions(${PLUGIN}::${NAME} PRIVATE
      # Deriving from std::iterator is deprecated in C++17 but used by VTK/ParaView. Revise when updating C++ standard!
      _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
      # std::unary_negate and std::binary_negate are deprecated in C++17 but used by Eigen. Revise when updating C++ standard!
      _SILENCE_CXX17_NEGATORS_DEPRECATION_WARNING
      # Many result_type typedefs are deprecated in C++17 but used by Eigen. Revise when updating C++ standard!
      _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
  endif()
endfunction()
