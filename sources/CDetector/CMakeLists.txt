add_subdirectory (CScoreProcessor)

add_library (detectors OBJECT CDetector.cpp)
target_link_libraries (detectors PUBLIC score_processor)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  set (LIST_DETECTOR_TB "")
  set (LIST_DETECTOR_DEP "")
  foreach (TARGET_Q ${LIST_TARGET_Q})
    set (LOCAL_FILENAME CDetector_q${TARGET_Q}_N60_w1_high_test.cpp)
    execute_process (COMMAND ${Python_EXECUTABLE} -c "import math as m; print(m.sqrt(${TARGET_Q}), end='')" OUTPUT_VARIABLE TARGET_SQRT_Q)
    list (APPEND LIST_DETECTOR_TB ${LOCAL_FILENAME})
    configure_file (testing/CDetector_qXX_N60_w1_high_test.cpp.in ${LOCAL_FILENAME} @ONLY)
    unset (LOCAL_FILENAME)
    unset (TARGET_SQRT_Q)

    list (APPEND LIST_DETECTOR_DEP gen_data_q${TARGET_Q}_N60_w1)
  endforeach ()

  add_executable (
    detector_tests testing/CDetector_w1_low_test.cpp testing/CDetector_w4_pi_4_low_test.cpp testing/CDetector_w34_pi_8_low_test.cpp
                   ${LIST_DETECTOR_TB}
  )
  target_link_libraries (detector_tests PUBLIC detectors catch_common_${PROJECT_NAME} ${matio})

  add_dependencies (detector_tests gen_data_q64_N60_w4 gen_data_q64_N60_w34 ${LIST_DETECTOR_DEP})

  catch_discover_tests (detector_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

endif ()
