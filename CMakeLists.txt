cmake_minimum_required (VERSION 3.16.0 FATAL_ERROR)
if (CMAKE_VERSION VERSION_LESS 3.18.0)
  message (WARNING "You are using an outdated version of cmake (typically used on Ubuntu 20.04). You may consider using an updated CMake.")
endif ()

# setting this is required

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set (CMAKE_EXPORT_COMPILE_COMMANDS true)

project (QCSP_PASSED LANGUAGES CXX VERSION 0.1)
# this sets the project name

if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set (IS_WINDOWS ON CACHE INTERNAL "System is Windows" FORCE)
  set (VCPKG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg CACHE INTERNAL "vcpkg directory" FORCE)
  message (AUTHOR_WARNING "Windows detected, adding vcpkg support.")
  set (ENV{VCPKG_DISABLE_METRICS} 1)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg/vcpkg.exe ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg/vcpkg.disable-metrics
    COMMAND ./bootstrap-vcpkg.bat -disableMetrics DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg/bootstrap-vcpkg.bat
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg
  )
  add_custom_target (
    vcpkg DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg/vcpkg.exe ${CMAKE_CURRENT_SOURCE_DIR}/windows/vcpkg/vcpkg.disable-metrics
  )

endif ()

add_compile_options (-Wall)
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  message (AUTHOR_WARNING "MSVC detected, adding specific options.")
  add_compile_options (/permissive-)
else ()
  add_compile_options (-Wextra)
endif ()

# ##########################################################################################################################################
# Options #####
# ##########################################################################################################################################

option (AGRESSIVE_OPTIMISATION "enable agressive flags" ON)
option (PREFER_STATIC_LIBRARY "prefer static libraries" ON)

option (PROFILING "enable GProf profiling" OFF)
option (DEBUG_SCORE "Dump score for debugging purpose" OFF)

option (ENABLE_TESTING "enable unit-testing with Catch2" ON)

include (CMakeDependentOption)

cmake_dependent_option (
  ENABLE_SYSTEM_CATCH2
  "Use system-provided Catch2 library"
  ON
  "NOT IS_WINDOWS;ENABLE_TESTING"
  OFF
)

cmake_dependent_option (
  ENABLE_DATA_AUTOGENERATION
  "Generate data automatically if an interpreter is found."
  OFF
  "ENABLE_TESTING"
  OFF
)

# ##########################################################################################################################################

include (FetchContent)

find_library (matio matio REQUIRED)
find_file (matio_HEADER matio.h REQUIRED)
get_filename_component (matio_INCLUDEDIR ${matio_HEADER} DIRECTORY CACHE)

# ## Currently unused dependencies

# find_package (Boost REQUIRED COMPONENTS program_options)

# find_library ( fftw NAMES fftw fftw3 HINTS "/usr/local/lib" "/usr/local" REQUIRED )

# find_library ( fftwf NAMES fftwf fftw3f HINTS "/usr/local/lib" "/usr/local" REQUIRED )

# find_package (OpenMP REQUIRED)

# find_package (Threads REQUIRED)

if (ENABLE_SYSTEM_CATCH2)
  find_package (Catch2 QUIET)
  if (((NOT Catch2_FOUND) OR (Catch2_VERSION VERSION_LESS 2.13)) AND ENABLE_TESTING)
    message (FATAL_ERROR "Catch2 not found, and you asked to use a system provided version. Please install Catch v2.x (x > 13)")
  endif ()
else ()
  fetchcontent_declare (
    Catch219 URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.9.tar.gz
    URL_HASH SHA256=06dbc7620e3b96c2b69d57bf337028bf245a211b3cddb843835bfe258f427a52
  )
  fetchcontent_populate (Catch219)
  add_subdirectory (${catch219_SOURCE_DIR} ${catch219_BINARY_DIR} EXCLUDE_FROM_ALL)
  list (APPEND CMAKE_MODULE_PATH "${catch219_SOURCE_DIR}/contrib")
endif ()

if (PROFILING)
  add_compile_definitions (WITHGPERFTOOLS)
  add_link_options (-lprofiler)
endif ()

include (CheckCXXCompilerFlag)
check_cxx_compiler_flag ("-march=native" MARCH)
if ((NOT MARCH))
  check_cxx_compiler_flag ("-mcpu=native" MCPU)
  if (MCPU)
    set (NATIVE_COMPILATION -mcpu=native)
  endif ()
else ()
  set (NATIVE_COMPILATION -march=native)
endif ()

check_cxx_compiler_flag ("-mtune=native" MTUNE)
if (MTUNE)
  set (NATIVE_COMPILATION ${NATIVE_COMPILATION} -mtune=native)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_definitions (DEBUG=2)
  add_compile_options (-O0 -ggdb ${NATIVE_COMPILATION})
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_definitions (NDEBUG)
  add_compile_options (-pedantic -Werror ${NATIVE_COMPILATION})
  if (AGRESSIVE_OPTIMISATION)
    add_compile_options (-Ofast)
  else ()
    add_compile_definitions (_FORTIFY_SOURCE=2)
  endif ()
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  add_compile_definitions (DEBUG=2)
  add_compile_options (-ggdb)
  add_compile_options (-pedantic -Werror ${NATIVE_COMPILATION})
  if (AGRESSIVE_OPTIMISATION)
    add_compile_options (-Ofast)
  else ()
    add_compile_definitions (_FORTIFY_SOURCE=2)
  endif ()
endif ()

if (DEBUG_SCORE)
  add_compile_definitions (DEBUG_LOG_SCORE=1)
endif ()

if (ENABLE_DATA_AUTOGENERATION)
  find_program (mat_interpreter NAMES matlab octave PATHS /bin /usr/bin /usr/local/bin QUIET)
  if (NOT mat_interpreter)
    message (STATUS "MATLAB/Octave interpreter not found. Test vectors will not be generated on-the-fly.")
  else ()
    string (REGEX MATCH "/matlab$" MAT_ITPR_TYPE ${mat_interpreter})
    if (MAT_ITPR_TYPE)
      if (NOT MATLAB_RELEASE)
        execute_process (COMMAND ${mat_interpreter} -nosplash -nojvm -batch "fprintf('%s', version())" OUTPUT_VARIABLE MATLAB_FUZZY)
        string (REGEX MATCH "R20[1-9][0-9][ab]" MATLAB_RELEASE_TMP ${MATLAB_FUZZY})
        set (MATLAB_RELEASE ${MATLAB_RELEASE_TMP} CACHE INTERNAL "" FORCE)
        unset (MATLAB_RELEASE_TMP)

        set (mat_generator_args "-nojvm;-nosplash;-batch" CACHE INTERNAL "")

        message (STATUS "Found MATLAB ${MATLAB_RELEASE} as ${mat_interpreter}.")
      endif ()

    else ()
      if (NOT OCTAVE_VERSION)
        execute_process (COMMAND ${mat_interpreter} --eval "fprintf('%s', OCTAVE_VERSION)" OUTPUT_VARIABLE OCTAVE_VERSION_TMP)

        set (OCTAVE_VERSION ${OCTAVE_VERSION_TMP} CACHE INTERNAL "" FORCE)
        unset (OCTAVE_VERSION_TMP)

        set (mat_generator_args "--no-gui;-f;-q;--eval" CACHE INTERNAL "")

        message (STATUS "Found Octave ${OCTAVE_VERSION} as  ${mat_interpreter}.")
      endif ()
    endif ()
  endif ()
endif ()

if (ENABLE_TESTING)
  find_package (Python COMPONENTS Interpreter REQUIRED)
endif ()

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/sources ${matio_INCLUDEDIR})

# ##########################################################################################################################################
# target definitions
# ##########################################################################################################################################
# ##########################################################################################################################################

set (LIST_TARGET_Q "64;128;256;512" CACHE STRING "Value of q to test.")
set (LIST_TARGET_N "60" CACHE INTERNAL "" FORCE) # Currently not supported
set (LIST_TARGET_W "1") # Currently not supported

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_library (catch_common_${PROJECT_NAME} OBJECT sources/testing/main_catch2.cpp)
  target_link_libraries (catch_common_${PROJECT_NAME} PUBLIC Catch2::Catch2)

  configure_file (data/generator_script.m.in ${CMAKE_CURRENT_SOURCE_DIR}/data/generator_script.m)

  foreach (TARGET_Q ${LIST_TARGET_Q})
    foreach (TARGET_N ${LIST_TARGET_N})
      foreach (TARGET_W ${LIST_TARGET_W})
        set (local_filename "test_data_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W}_step0_span0.0.mat")

        configure_file (data/data_generator_qXX_NXX_wXX_step0_0pi_1.m.in ${CMAKE_CURRENT_SOURCE_DIR}/data/data_generator_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W}_step0_0pi_1.m)

        # If MATLAB or Octave have been found, generate test vectors. Just create dummy dependencies otherwise.
        if ((OCTAVE_VERSIONS OR MATLAB_RELEASES) AND ENABLE_DATA_AUTOGENERATION)
          add_custom_command (
            OUTPUT ${CMAKE_SOURCE_DIR}/data/${local_filename}
            COMMAND
              ${mat_interpreter} ${mat_generator_args}
              "PN = double(load('pn_sequences.mat').PN${TARGET_Q}); best_N = load('parameters_20210903.mat').best_N; fnm = save_test_vectors(PN, best_N, 0, 1, 0, 0); fprintf('%s generated.\\n', fnm)"
            DEPENDS ${CMAKE_SOURCE_DIR}/data/save_test_vectors.m
                    ${CMAKE_SOURCE_DIR}/data/compute_ts_score.m
                    ${CMAKE_SOURCE_DIR}/data/generate_noisy_sequence.m
                    ${CMAKE_SOURCE_DIR}/data/corr_abs_max.m
                    ${CMAKE_SOURCE_DIR}/data/pn_sequences.mat
                    ${CMAKE_SOURCE_DIR}/data/parameters_20210903.mat
            USES_TERMINAL
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data
            VERBATIM
          )
        endif ()

        add_custom_target (gen_data_q${TARGET_Q}_N${TARGET_N}_w${TARGET_W} DEPENDS ${CMAKE_SOURCE_DIR}/data/${local_filename})
        unset (${local_filename})
      endforeach ()
    endforeach ()
  endforeach ()

endif ()

add_subdirectory (sources/CDetector)
