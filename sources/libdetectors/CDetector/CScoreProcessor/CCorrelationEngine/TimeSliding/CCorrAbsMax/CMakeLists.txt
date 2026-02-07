add_library (corrabsmax STATIC CCorrAbsMax.cpp)

add_library (corrabsmax_generic STATIC CCorrAbsMaxGeneric.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (corrabsmax_tests testing/CCorrAbsMax_test.cpp testing/CCorrAbsMax_i16_test.cpp)
  target_link_libraries (corrabsmax_tests PUBLIC corrabsmax catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (corrabsmax_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (corrabsmax_generic_tests testing/CCorrAbsMaxGeneric_test.cpp)
  target_link_libraries (corrabsmax_generic_tests PUBLIC corrabsmax_generic catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (corrabsmax_generic_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
