add_library (corrabsmax OBJECT CCorrAbsMax.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (corrabsmax_tests testing/CCorrAbsMax_test.cpp testing/CCorrAbsMax_i16_test.cpp)
  target_link_libraries (corrabsmax_tests PUBLIC corrabsmax catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (corrabsmax_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
