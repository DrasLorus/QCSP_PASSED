add_library (qspannedsequentialadder STATIC CQSpannedSequentialAdder.cpp)

add_library (qspannedsequentialadder_rc STATIC CQSpannedSequentialAdderRC.cpp)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  add_executable (qspannedsequentialadder_tests testing/CQSpannedSequentialAdder_test.cpp testing/CQSpannedSequentialAdder_i16_test.cpp)
  target_link_libraries (qspannedsequentialadder_tests PUBLIC qspannedsequentialadder catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (qspannedsequentialadder_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable (qspannedsequentialadder_rc_tests testing/CQSpannedSequentialAdderRC_test.cpp)
  target_link_libraries (qspannedsequentialadder_rc_tests PUBLIC qspannedsequentialadder_rc catch_common_${PROJECT_NAME} ${matio})

  catch_discover_tests (qspannedsequentialadder_rc_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
