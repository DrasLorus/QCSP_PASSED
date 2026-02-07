cmake_minimum_required (VERSION 3.16.0 FATAL_ERROR)

project (VCPKG_WRAPPER VERSION 0.1)

set (IS_WINDOWS ON CACHE INTERNAL "System is Windows" FORCE)
set (VCPKG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg CACHE INTERNAL "vcpkg directory" FORCE)
set (ENV{VCPKG_DISABLE_METRICS} 1)

add_custom_command (
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/vcpkg.exe ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/vcpkg.disable-metrics
  COMMAND ./bootstrap-vcpkg.bat -disableMetrics
  DEPENDS ${CMAKE_SOURCE_DIR}/vcpkg.json ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/bootstrap-vcpkg.bat
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg
)

add_custom_target (
  vcpkg DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/vcpkg.exe ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/vcpkg.disable-metrics
                ${CMAKE_SOURCE_DIR}/vcpkg.json
)
