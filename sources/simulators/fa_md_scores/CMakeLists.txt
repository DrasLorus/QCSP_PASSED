find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS program_options REQUIRED)

include (CheckIncludeFiles)
check_include_files ("unistd.h" HAVE_UNISTD_H)

configure_file (config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_executable (fa_md_simulator main.cpp)
target_link_libraries (fa_md_simulator PRIVATE Boost::program_options)
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_link_libraries (fa_md_simulator PRIVATE Pathcch.lib)
endif ()

set_target_properties (fa_md_simulator PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_library (fa_md_utilities OBJECT utility_functions.cpp)
target_link_libraries (fa_md_utilities PUBLIC ${matio})
target_include_directories (fa_md_utilities PUBLIC ${matio_INCLUDEDIR})

set (FAMD_TARGET_N "60")
set (FAMD_TARGET_Q "64;256")
set (FAMD_TARGET_W "1;2;3;4;5;8")
set (FAMD_TARGET_NORMED "0;1")

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set (FAMD_LOCATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/fa_md_simulation)
else ()
  set (FAMD_LOCATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fa_md_simulation)
endif ()

unset (USE_MULT)
foreach (TARGET_N ${FAMD_TARGET_N})
  foreach (TARGET_Q ${FAMD_TARGET_Q})
    foreach (TARGET_W ${FAMD_TARGET_W})
      foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})
        set (LOCAL_SUFFIX N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED})

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
        )

        set_target_properties (fa_md_simulation_${LOCAL_SUFFIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${FAMD_LOCATION})

        add_dependencies (fa_md_simulator fa_md_simulation_${LOCAL_SUFFIX})
      endforeach ()
    endforeach ()
  endforeach ()
endforeach ()

foreach (TARGET_N ${FAMD_TARGET_N})
  foreach (TARGET_Q ${FAMD_TARGET_Q})
    foreach (TARGET_W ${FAMD_TARGET_W})
      foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})
        set (USE_MULT 1)
        set (LOCAL_SUFFIX N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}_mult)

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
        )

        set_target_properties (fa_md_simulation_${LOCAL_SUFFIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${FAMD_LOCATION})

        add_dependencies (fa_md_simulator fa_md_simulation_${LOCAL_SUFFIX})
      endforeach ()
    endforeach ()
  endforeach ()
endforeach ()
