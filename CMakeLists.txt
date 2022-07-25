cmake_minimum_required (VERSION 3.18.0 FATAL_ERROR)
# setting this is required

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set (CMAKE_EXPORT_COMPILE_COMMANDS true)

project (
  QCSP_PASSED
  LANGUAGES CXX
  VERSION 0.1
)
# this sets the project name

# ##################################################################################################
# Options #####
# ##################################################################################################

option (AGRESSIVE_OPTIMISATION "enable agressive flags" ON)
option (PREFER_STATIC_LIBRARY "prefer static libraries" ON)

option (PROFILING "enable GProf profiling" OFF)
option (DEBUG_SCORE "Dump score for debugging purpose" OFF)

option (ENABLE_TESTING "enable unit-testing with Catch2" ON)

# ##################################################################################################

find_package (Boost REQUIRED COMPONENTS program_options)

find_library (matio matio REQUIRED)

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

add_library (l2norm OBJECT sources/CDetector/CScoreProcessor/CNorm/CNorm.cpp)
target_include_directories (l2norm PUBLIC sources/CDetector sources)

add_library (
  iterativeadder OBJECT sources/CDetector/CScoreProcessor/CIterativeAdder/CIterativeAdder.cpp
)
target_include_directories (iterativeadder PUBLIC sources/CDetector sources)

add_library (corrabsmax OBJECT sources/CDetector/CScoreProcessor/CCorrAbsMax/CCorrAbsMax.cpp)
target_include_directories (corrabsmax PUBLIC sources/CDetector sources)

add_library (
  scoreaccu OBJECT sources/CDetector/CScoreProcessor/CScoreAccumulator/CScoreAccumulator.cpp
)
target_include_directories (scoreaccu PUBLIC sources/CDetector sources)

add_library (score_processor OBJECT sources/CDetector/CScoreProcessor/CScoreProcessor.cpp)
target_include_directories (score_processor PUBLIC sources sources/CDetector)
target_link_libraries (
  score_processor
  PUBLIC corrabsmax
         l2norm
         iterativeadder
         scoreaccu
)

add_library (detectors OBJECT sources/CDetector/CDetector.cpp)
target_include_directories (score_processor PUBLIC sources)
target_link_libraries (detectors PUBLIC score_processor)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_library (catch_common_${PROJECT_NAME} OBJECT sources/testing/main_catch2.cpp)
  target_link_libraries (catch_common_${PROJECT_NAME} PUBLIC Catch2::Catch2)

  add_executable (
    l2norm_tests sources/CDetector/CScoreProcessor/CNorm/testing/CNorm_test.cpp
                 sources/CDetector/CScoreProcessor/CNorm/testing/CNorm_i16_test.cpp
  )
  target_link_libraries (l2norm_tests PUBLIC l2norm catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (l2norm_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    iterativeadder_tests
    sources/CDetector/CScoreProcessor/CIterativeAdder/testing/CIterativeAdder_test.cpp
    sources/CDetector/CScoreProcessor/CIterativeAdder/testing/CIterativeAdder_i16_test.cpp
  )
  target_link_libraries (
    iterativeadder_tests PUBLIC iterativeadder catch_common_${PROJECT_NAME} matio
  )

  catch_discover_tests (iterativeadder_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    corrabsmax_tests sources/CDetector/CScoreProcessor/CCorrAbsMax/testing/CCorrAbsMax_test.cpp
                     sources/CDetector/CScoreProcessor/CCorrAbsMax/testing/CCorrAbsMax_i16_test.cpp
  )
  target_link_libraries (corrabsmax_tests PUBLIC corrabsmax catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (corrabsmax_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    scoreaccu_tests
    sources/CDetector/CScoreProcessor/CScoreAccumulator/testing/CScoreAccumulator_test.cpp
    sources/CDetector/CScoreProcessor/CScoreAccumulator/testing/CScoreAccumulator_i16_test.cpp
  )
  target_link_libraries (scoreaccu_tests PUBLIC scoreaccu catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (scoreaccu_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    scoreproc_tests sources/CDetector/CScoreProcessor/testing/CScoreProcessor_test.cpp
                    sources/CDetector/CScoreProcessor/testing/CScoreProcessor_i16_test.cpp
  )
  target_link_libraries (scoreproc_tests PUBLIC score_processor catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (scoreproc_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (
    detector_tests sources/CDetector/testing/CDetector_w1_high_test.cpp
                   sources/CDetector/testing/CDetector_w1_low_test.cpp
  )
  target_link_libraries (detector_tests PUBLIC detectors catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (detector_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
