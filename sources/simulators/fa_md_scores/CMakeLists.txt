find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS program_options REQUIRED)

include (CheckIncludeFiles)
check_include_files ("unistd.h" HAVE_UNISTD_H)
check_include_files ("windows.h" HAVE_WINDOWS_H)

if (WIN32 AND HAVE_WINDOWS_H)
  set (USE_UNISTD_API OFF)
  set (USE_WINDOWS_API ON)
elseif (HAVE_UNISTD_H)
  set (USE_UNISTD_API ON)
  set (USE_WINDOWS_API OFF)
endif ()

option (ENABLE_FAMD_MULT "Wether or not you want to instantiate plain multipliers. No benefit." OFF)
set (USE_MULT ${ENABLE_FAMD_MULT} CACHE INTERNAL "" FORCE)

configure_file (config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_executable (fa_md_simulator main.cpp)
target_link_libraries (fa_md_simulator PRIVATE Boost::program_options)
if (USE_WINDOWS_API)
  target_link_libraries (fa_md_simulator PRIVATE Pathcch.lib)
endif ()

set_target_properties (fa_md_simulator PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_library (fa_md_utilities OBJECT utility_functions.cpp)
target_link_libraries (fa_md_utilities PUBLIC ${matio})
target_include_directories (fa_md_utilities PUBLIC ${matio_INCLUDEDIR})

set (FAMD_TARGET_N "60" CACHE INTERNAL "Not supported yet")
set (FAMD_TARGET_Q "64;256" CACHE STRING "List of GF(q) to support. Warning: corresponding AList and PN must be available.")
set (FAMD_TARGET_W "0;1;2;3;4;5;8" CACHE STRING "List of p_omega to support.")
set (FAMD_TARGET_NORMED "0;1" CACHE STRING "List of norm method to support. Either normed (1) or not (0).")
# set (FAMD_TARGET_INTYPE "float;int16_t" CACHE STRING "List of type to support. Only 'float' and 'int16_t' are supported, using something
# else is undefined behavior." )

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set (FAMD_LOCATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/fa_md_simulation)
else ()
  set (FAMD_LOCATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fa_md_simulation)
endif ()

if (NOT USE_MULT)
  set (TARGET_INTYPE "float")
  foreach (TARGET_N ${FAMD_TARGET_N})
    foreach (TARGET_Q ${FAMD_TARGET_Q})
      foreach (TARGET_W ${FAMD_TARGET_W})
        foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})
          set (LOCAL_SUFFIX N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_INTYPE}_${TARGET_NORMED})

          configure_file (detector_instantiation.cpp.in detector_${LOCAL_SUFFIX}.cpp @ONLY)

          add_library (libdetector_${LOCAL_SUFFIX} OBJECT detector_${LOCAL_SUFFIX}.cpp)
          target_link_libraries (libdetector_${LOCAL_SUFFIX} PRIVATE detectors)

          configure_file (main_N_q_w_normed.cpp.in main_${LOCAL_SUFFIX}.cpp @ONLY)

          add_executable (fa_md_simulation_${LOCAL_SUFFIX} main_${LOCAL_SUFFIX}.cpp)

          target_include_directories (fa_md_simulation_${LOCAL_SUFFIX} BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
          target_link_directories (fa_md_simulation_${LOCAL_SUFFIX} PRIVATE ${matio_INCLUDEDIR})
          target_link_libraries (
            fa_md_simulation_${LOCAL_SUFFIX}
            PRIVATE generator
                    detectors
                    fa_md_utilities
                    Boost::program_options
                    OpenMP::OpenMP_CXX
                    ${matio}
                    nbldpc_wrap
                    libdetector_${LOCAL_SUFFIX}
          )

          set_target_properties (fa_md_simulation_${LOCAL_SUFFIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${FAMD_LOCATION})

          add_dependencies (fa_md_simulator fa_md_simulation_${LOCAL_SUFFIX})
        endforeach ()
      endforeach ()
    endforeach ()
  endforeach ()

  set (TARGET_INTYPE "int16_t")
  foreach (TARGET_N ${FAMD_TARGET_N})
    foreach (TARGET_Q ${FAMD_TARGET_Q})
      foreach (TARGET_W ${FAMD_TARGET_W})
        foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})
          set (LOCAL_SUFFIX N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_INTYPE}_${TARGET_NORMED})

          configure_file (detector_instantiation_int16.cpp.in detector_${LOCAL_SUFFIX}.cpp @ONLY)

          add_library (libdetector_${LOCAL_SUFFIX} OBJECT detector_${LOCAL_SUFFIX}.cpp)
          target_link_libraries (libdetector_${LOCAL_SUFFIX} PRIVATE detectors)

          configure_file (main_N_q_w_normed_int16.cpp.in main_${LOCAL_SUFFIX}.cpp @ONLY)

          add_executable (fa_md_simulation_${LOCAL_SUFFIX} main_${LOCAL_SUFFIX}.cpp)

          target_include_directories (fa_md_simulation_${LOCAL_SUFFIX} BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
          target_link_directories (fa_md_simulation_${LOCAL_SUFFIX} PRIVATE ${matio_INCLUDEDIR})
          target_link_libraries (
            fa_md_simulation_${LOCAL_SUFFIX}
            PRIVATE generator
                    detectors
                    fa_md_utilities
                    Boost::program_options
                    OpenMP::OpenMP_CXX
                    ${matio}
                    nbldpc_wrap
                    libdetector_${LOCAL_SUFFIX}
          )

          set_target_properties (fa_md_simulation_${LOCAL_SUFFIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${FAMD_LOCATION})

          add_dependencies (fa_md_simulator fa_md_simulation_${LOCAL_SUFFIX})
        endforeach ()
      endforeach ()
    endforeach ()
  endforeach ()
else ()
  set (TARGET_INTYPE "float")
  foreach (TARGET_N ${FAMD_TARGET_N})
    foreach (TARGET_Q ${FAMD_TARGET_Q})
      foreach (TARGET_W ${FAMD_TARGET_W})
        foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})
          set (LOCAL_SUFFIX N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_INTYPE}_${TARGET_NORMED}_mult)

          configure_file (detector_instantiation.cpp.in detector_${LOCAL_SUFFIX}.cpp @ONLY)

          add_library (libdetector_${LOCAL_SUFFIX} OBJECT detector_${LOCAL_SUFFIX}.cpp)
          target_link_libraries (libdetector_${LOCAL_SUFFIX} PRIVATE detectors)

          configure_file (main_N_q_w_normed.cpp.in main_${LOCAL_SUFFIX}.cpp @ONLY)

          add_executable (fa_md_simulation_${LOCAL_SUFFIX} main_${LOCAL_SUFFIX}.cpp)

          target_include_directories (fa_md_simulation_${LOCAL_SUFFIX} BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
          target_link_directories (fa_md_simulation_${LOCAL_SUFFIX} PRIVATE ${matio_INCLUDEDIR})
          target_link_libraries (
            fa_md_simulation_${LOCAL_SUFFIX}
            PRIVATE generator
                    detectors
                    fa_md_utilities
                    Boost::program_options
                    OpenMP::OpenMP_CXX
                    ${matio}
                    nbldpc_wrap
                    libdetector_${LOCAL_SUFFIX}
          )

          set_target_properties (fa_md_simulation_${LOCAL_SUFFIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${FAMD_LOCATION})

          add_dependencies (fa_md_simulator fa_md_simulation_${LOCAL_SUFFIX})
        endforeach ()
      endforeach ()
    endforeach ()
  endforeach ()
endif ()

unset (USE_MULT)
