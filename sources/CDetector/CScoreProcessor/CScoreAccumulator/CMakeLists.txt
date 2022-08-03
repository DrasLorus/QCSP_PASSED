add_library (scoreaccu OBJECT CScoreAccumulator.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (
    scoreaccu_tests testing/CScoreAccumulator_test.cpp testing/CScoreAccumulator_i16_test.cpp
  )
  target_link_libraries (scoreaccu_tests PUBLIC scoreaccu catch_common_${PROJECT_NAME} matio)

  catch_discover_tests (scoreaccu_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
