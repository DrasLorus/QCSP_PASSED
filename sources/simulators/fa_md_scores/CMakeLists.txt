find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS program_options REQUIRED)

set (FAMD_TARGET_N "60")
set (FAMD_TARGET_Q "64")
set (FAMD_TARGET_W "1;2;4;8")
set (FAMD_TARGET_NORMED "0;1")

foreach (TARGET_N ${FAMD_TARGET_N})
  foreach (TARGET_Q ${FAMD_TARGET_Q})
    foreach (TARGET_W ${FAMD_TARGET_W})
      foreach (TARGET_NORMED ${FAMD_TARGET_NORMED})

        configure_file (main_N_q_w_normed.cpp.in main_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}.cpp @ONLY)

        add_executable (
          fa_md_simulation_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}
          main_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}.cpp
        )

        target_link_directories (fa_md_simulation_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED} PRIVATE ${matio_INCLUDEDIR})
        target_link_libraries (
          fa_md_simulation_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}
          PRIVATE generator
                  detectors
                  Boost::program_options
                  OpenMP::OpenMP_CXX
                  ${matio}
                  nbldpc_wrap
        )

        set_target_properties (
          fa_md_simulation_N${TARGET_N}_q${TARGET_Q}_w${TARGET_W}_${TARGET_NORMED}
          PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fa_md_simulation
        )
      endforeach ()
    endforeach ()
  endforeach ()
endforeach ()

include (CheckIncludeFiles)
check_include_files ("unistd.h" HAVE_UNISTD_H)

configure_file (config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_executable (fa_md_simulator main.cpp)
target_link_libraries (fa_md_simulator PRIVATE Boost::program_options)

set_target_properties (fa_md_simulator PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
