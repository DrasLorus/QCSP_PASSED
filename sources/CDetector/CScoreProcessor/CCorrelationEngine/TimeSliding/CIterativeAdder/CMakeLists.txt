add_library (iterativeadder OBJECT CIterativeAdder.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (iterativeadder_tests testing/CIterativeAdder_test.cpp testing/CIterativeAdder_i16_test.cpp)
  target_link_libraries (iterativeadder_tests PUBLIC iterativeadder catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (iterativeadder_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
