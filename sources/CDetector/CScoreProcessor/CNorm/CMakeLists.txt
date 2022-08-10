add_library (l2norm OBJECT CNorm.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (l2norm_tests testing/CNorm_test.cpp testing/CNorm_i16_test.cpp)
  target_link_libraries (l2norm_tests PUBLIC l2norm catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (l2norm_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
