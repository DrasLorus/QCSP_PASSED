add_subdirectory (CNorm)

add_subdirectory (CCorrelationEngine)

add_subdirectory (CScoreAccumulator)

add_library (score_processor STATIC CScoreProcessor.cpp)
target_link_libraries (score_processor PUBLIC l2norm scoreaccu correngine)

add_library (score_processor_generic STATIC CScoreProcessorGeneric.cpp)
target_link_libraries (score_processor_generic PUBLIC l2norm_generic scoreaccu_generic correngine_generic)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (scoreproc_tests testing/CScoreProcessor_test.cpp testing/CScoreProcessor_i16_test.cpp)
  target_link_libraries (scoreproc_tests PUBLIC score_processor catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (scoreproc_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (scoreproc_generic_tests testing/CScoreProcessorGeneric_test.cpp CScoreProcessorGeneric.cpp)
  target_link_libraries (scoreproc_generic_tests PUBLIC score_processor_generic catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (scoreproc_generic_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
