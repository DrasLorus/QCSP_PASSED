cmake_minimum_required (VERSION 3.18.0 FATAL_ERROR)
# setting this is required

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set (CMAKE_EXPORT_COMPILE_COMMANDS true)

project (
  QCSP_LoneDetector
  LANGUAGES CXX
  VERSION 0.2
)
# this sets the project name

# ##################################################################################################
# Options #####
# ##################################################################################################

option (AGRESSIVE_OPTIMISATION "enable agressive flags" ON)
option (PREFER_STATIC_LIBRARY "prefer static libraries" ON)

option (PROFILING "enable GProf profiling" OFF)
option (DEBUG_SCORE "Dump score for debugging purpose" OFF)

option (ENABLE_TESTING "enable unit-testing with Catch2" OFF)

# ##################################################################################################

if (PREFER_STATIC_LIBRARY)
  set (BLA_STATIC ON)
endif ()

find_package (BLAS REQUIRED)
message (VERBOSE "BLAS Vendor: ${BLA_VENDOR}. Found: ${BLAS_LIBRARIES}")

find_package (Boost REQUIRED COMPONENTS program_options)

find_package (matioCpp QUIET)
if (matioCpp_FOUND EQUAL 0)
  set (matioCpp_VERSION 0.2.0)
  message (STATUS "matioCpp not found. Version ${matioCpp_VERSION} will be build locally.")
  file (
    DOWNLOAD "https://github.com/dic-iit/matio-cpp/archive/refs/tags/v${matioCpp_VERSION}.zip"
    "${CMAKE_BINARY_DIR}/matio-cpp-${matioCpp_VERSION}.zip"
    EXPECTED_HASH SHA256=c102638834e14b923053fb4f0baf3f258ad2a055dc2bc88e658725c4a72789e1
    TLS_VERIFY ON
  )
  if (CMAKE_VERSION VERSION_LESS 3.18)
    message (STATUS "Please consider to switch to CMake 3.18+")
    execute_process (
      COMMAND "${CMAKE_COMMAND}" -E tar xvzf "${CMAKE_BINARY_DIR}/matio-cpp-${matioCpp_VERSION}.zip"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
  else ()
    file (
      ARCHIVE_EXTRACT
      INPUT
      "${CMAKE_BINARY_DIR}/matio-cpp-${matioCpp_VERSION}.zip"
      DESTINATION
      "${CMAKE_BINARY_DIR}"
    )
  endif ()
  add_subdirectory ("${CMAKE_BINARY_DIR}/matio-cpp-${matioCpp_VERSION}")
else ()
  message (STATUS "Found matioCpp: ${matioCpp_INCLUDE_DIRS}")
endif ()

find_library (
  fftw
  NAMES fftw fftw3
  HINTS "/usr/local/lib" "/usr/local" REQUIRED
)

find_library (
  fftwf
  NAMES fftwf fftw3f
  HINTS "/usr/local/lib" "/usr/local" REQUIRED
)

find_package (Catch2 QUIET)
if (NOT Catch2_FOUND)
  message (WARNING "Catch2 not found, unit tests and checking are disabled.")
endif ()

find_package (OpenMP REQUIRED)
if (OPENMP_FFT_GRID)
  add_compile_definitions (OMP_FFT_P)
endif ()

find_package (Threads REQUIRED)

if (PROFILING)
  add_compile_definitions (WITHGPERFTOOLS)
  add_link_options (-lprofiler)
endif ()

add_compile_options (-Wall -Wextra)

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
  add_compile_options (-O0 -ggdb)
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

# ##################################################################################################
# target definitions
# ##################################################################################################
# ##################################################################################################

add_library (l2norm OBJECT sources/CScoreProcessor/CNorm/CNorm.cpp)
target_include_directories (l2norm PUBLIC sources)

add_library (iterativeadder OBJECT sources/CScoreProcessor/CIterativeAdder/CIterativeAdder.cpp)
target_include_directories (iterativeadder PUBLIC sources)

add_library (corrabsmax OBJECT sources/CScoreProcessor/CCorrAbsMax/CCorrAbsMax.cpp)
target_include_directories (corrabsmax PUBLIC sources)

add_library (scoreaccu OBJECT sources/CScoreProcessor/CScoreAccumulator/CScoreAccumulator.cpp)
target_include_directories (scoreaccu PUBLIC sources)

add_library (score_processor OBJECT sources/CScoreProcessor/CScoreProcessor.cpp)
target_include_directories (score_processor PUBLIC sources sources/CScoreProcessor)
target_link_libraries (
  score_processor
  PUBLIC corrabsmax
         l2norm
         iterativeadder
         scoreaccu
)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_library (catch_common_${PROJECT_NAME} OBJECT sources/testing/main_catch2.cpp)
  target_link_libraries (catch_common_${PROJECT_NAME} PUBLIC Catch2::Catch2)

  add_executable (
    l2norm_tests sources/CScoreProcessor/CNorm/testing/CNorm_test.cpp
                 sources/CScoreProcessor/CNorm/testing/CNorm_i16_test.cpp
  )
  target_link_libraries (l2norm_tests PUBLIC l2norm catch_common_${PROJECT_NAME} matioCpp::matioCpp)

  catch_discover_tests (l2norm_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    iterativeadder_tests
    sources/CScoreProcessor/CIterativeAdder/testing/CIterativeAdder_test.cpp
    sources/CScoreProcessor/CIterativeAdder/testing/CIterativeAdder_i16_test.cpp
  )
  target_link_libraries (
    iterativeadder_tests PUBLIC iterativeadder catch_common_${PROJECT_NAME} matioCpp::matioCpp
  )

  catch_discover_tests (iterativeadder_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    corrabsmax_tests sources/CScoreProcessor/CCorrAbsMax/testing/CCorrAbsMax_test.cpp
                     sources/CScoreProcessor/CCorrAbsMax/testing/CCorrAbsMax_i16_test.cpp
  )
  target_link_libraries (
    corrabsmax_tests PUBLIC corrabsmax catch_common_${PROJECT_NAME} matioCpp::matioCpp
  )

  catch_discover_tests (corrabsmax_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    scoreaccu_tests
    sources/CScoreProcessor/CScoreAccumulator/testing/CScoreAccumulator_test.cpp
    sources/CScoreProcessor/CScoreAccumulator/testing/CScoreAccumulator_i16_test.cpp
  )
  target_link_libraries (
    scoreaccu_tests PUBLIC scoreaccu catch_common_${PROJECT_NAME} matioCpp::matioCpp
  )

  catch_discover_tests (scoreaccu_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    scoreproc_tests sources/CScoreProcessor/testing/CScoreProcessor_test.cpp
                    sources/CScoreProcessor/testing/CScoreProcessor_i16_test.cpp
  )
  target_link_libraries (
    scoreproc_tests PUBLIC scoreaccu catch_common_${PROJECT_NAME} matioCpp::matioCpp
  )

  catch_discover_tests (scoreproc_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
