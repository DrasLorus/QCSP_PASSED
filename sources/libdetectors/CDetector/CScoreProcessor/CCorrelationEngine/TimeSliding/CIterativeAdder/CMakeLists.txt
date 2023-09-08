add_library (iterativeadder STATIC CIterativeAdder.cpp)

add_library (iterativeadder_generic STATIC CIterativeAdderGeneric.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (iterativeadder_tests testing/CIterativeAdder_test.cpp testing/CIterativeAdder_i16_test.cpp)
  target_link_libraries (iterativeadder_tests PUBLIC iterativeadder catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (iterativeadder_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (iterativeadder_generic_tests testing/CIterativeAdderGeneric_test.cpp)
  target_link_libraries (iterativeadder_generic_tests PUBLIC iterativeadder_generic catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (iterativeadder_generic_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
