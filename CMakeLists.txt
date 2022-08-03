cmake_minimum_required (VERSION 3.16.0 FATAL_ERROR)
if (CMAKE_VERSION VERSION_LESS 3.18.0)
  message (
    WARNING
      "You are using an outdated version of cmake (typically used on Ubuntu 20.04). You may consider using an updated CMake."
  )
endif ()

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

option (ENABLE_SYSTEM_CATCH2 "Use system-provided Catch2 library" ON)

# ##################################################################################################

include (FetchContent)

find_package (Boost REQUIRED COMPONENTS program_options)

find_library (matio matio REQUIRED)
find_file (matio-header matio.h REQUIRED)

# Currently unused dependencies

# find_library ( fftw NAMES fftw fftw3 HINTS "/usr/local/lib" "/usr/local" REQUIRED )

# find_library ( fftwf NAMES fftwf fftw3f HINTS "/usr/local/lib" "/usr/local" REQUIRED )

# find_package (OpenMP REQUIRED)

# find_package (Threads REQUIRED)

find_package (Catch2 QUIET)
if (((NOT Catch2_FOUND) OR (Catch2_VERSION VERSION_LESS 2.13)) AND ENABLE_TESTING)
  if (ENABLE_SYSTEM_CATCH2)
    message (
      FATAL_ERROR
        "Catch2 not found, and you asked to use a system provided version. Please install Catch v2.x (x > 13)"
    )
  else ()
    fetchcontent_declare (
      Catch219
      URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.9.tar.gz
      URL_HASH SHA256=06dbc7620e3b96c2b69d57bf337028bf245a211b3cddb843835bfe258f427a52
    )
    fetchcontent_makeavailable (Catch219)
    list (APPEND CMAKE_MODULE_PATH "${catch219_SOURCE_DIR}/contrib")
  endif ()
endif ()

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

include_directories (BEFORE sources)

# ##################################################################################################
# target definitions
# ##################################################################################################
# ##################################################################################################

set (LIST_TARGET_Q "64;128;256;512")

add_subdirectory (sources/CDetector)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_library (catch_common_${PROJECT_NAME} OBJECT sources/testing/main_catch2.cpp)
  target_link_libraries (catch_common_${PROJECT_NAME} PUBLIC Catch2::Catch2)

endif ()
