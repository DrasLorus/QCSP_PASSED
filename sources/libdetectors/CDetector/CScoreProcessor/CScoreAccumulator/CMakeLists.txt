add_library (scoreaccu STATIC CScoreAccumulator.cpp)

add_library (scoreaccu_generic STATIC CScoreAccumulatorGeneric.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (scoreaccu_tests testing/CScoreAccumulator_test.cpp testing/CScoreAccumulator_i16_test.cpp)
  target_link_libraries (scoreaccu_tests PUBLIC scoreaccu catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (scoreaccu_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (scoreaccu_generic_tests testing/CScoreAccumulatorGeneric_test.cpp)
  target_link_libraries (scoreaccu_generic_tests PUBLIC scoreaccu_generic catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (scoreaccu_generic_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
