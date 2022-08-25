add_subdirectory (CNorm)

add_subdirectory (CCorrelationEngine)

add_subdirectory (CScoreAccumulator)

add_library (score_processor OBJECT CScoreProcessor.cpp)

target_link_libraries (score_processor PUBLIC l2norm scoreaccu correngine)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (scoreproc_tests testing/CScoreProcessor_test.cpp testing/CScoreProcessor_i16_test.cpp)
  target_link_libraries (scoreproc_tests PUBLIC score_processor catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (scoreproc_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
