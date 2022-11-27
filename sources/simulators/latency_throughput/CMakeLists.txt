find_package (Threads REQUIRED)
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

configure_file (config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

foreach (TARGET_N ${LIST_TARGET_N})
  foreach (TARGET_Q ${LIST_TARGET_Q})
    foreach (TARGET_W ${LIST_TARGET_W})
      configure_file (main_benchmark.cpp.in detector_benchmark_${TARGET_N}_${TARGET_Q}_${TARGET_W}.cpp @ONLY)

      add_executable (detector_benchmark_${TARGET_N}_${TARGET_Q}_${TARGET_W} detector_benchmark_${TARGET_N}_${TARGET_Q}_${TARGET_W}.cpp)

      target_include_directories (detector_benchmark_${TARGET_N}_${TARGET_Q}_${TARGET_W} BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
      target_link_libraries (
        detector_benchmark_${TARGET_N}_${TARGET_Q}_${TARGET_W} PRIVATE detectors ${matio} Boost::program_options Threads::Threads
      )
    endforeach ()
  endforeach ()
endforeach ()

add_executable (detector_benchmark main.cpp)

target_include_directories (detector_benchmark BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (detector_benchmark PRIVATE detectors ${matio} Boost::program_options Threads::Threads)
