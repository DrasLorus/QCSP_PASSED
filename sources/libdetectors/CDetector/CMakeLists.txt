add_subdirectory (CScoreProcessor)

add_library (detectors STATIC CDetector.cpp)
target_include_directories (detectors PUBLIC ${LIB_ROOT_DIRECTORY})
target_link_libraries (detectors PUBLIC score_processor)

add_library (detectors_generic STATIC CDetectorSerialGeneric.cpp CDetectionStateGeneric.cpp)
target_include_directories (detectors_generic PUBLIC ${LIB_ROOT_DIRECTORY})
target_link_libraries (detectors_generic PUBLIC score_processor_generic)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  # For templated detectors
  set (LIST_DETECTOR_TB "")
  set (LIST_DETECTOR_DEP "")
  foreach (TARGET_Q ${LIST_TARGET_Q})
    execute_process (COMMAND ${Python_EXECUTABLE} -c "import math as m; print(m.sqrt(${TARGET_Q}), end='')" OUTPUT_VARIABLE TARGET_SQRT_Q)

    set (LOCAL_FILENAME CDetector_q${TARGET_Q}_N60_w1_high_test.cpp)
    list (APPEND LIST_DETECTOR_TB ${LOCAL_FILENAME})
    configure_file (testing/CDetector_qXX_N60_w1_high_test.cpp.in ${LOCAL_FILENAME} @ONLY)
    unset (LOCAL_FILENAME)

    set (LOCAL_FILENAME CDetectorSerialMult_q${TARGET_Q}_N60_w1_high_test.cpp)
    list (APPEND LIST_DETECTOR_TB ${LOCAL_FILENAME})
    configure_file (testing/CDetectorSerialMult_qXX_N60_w1_high_test.cpp.in ${LOCAL_FILENAME} @ONLY)
    unset (LOCAL_FILENAME)

    list (APPEND LIST_DETECTOR_DEP gen_data_q${TARGET_Q}_N60_w1)
    unset (TARGET_SQRT_Q)
  endforeach ()

  add_executable (
    detectors_tests
    testing/CDetector_w1_low_test.cpp
    testing/CDetector_w4_pi_4_low_test.cpp
    testing/CDetector_w34_pi_8_low_test.cpp
    testing/CDetectorSerialMult_w34_pi_8_low_test.cpp
    ${LIST_DETECTOR_TB}
  )
  target_link_libraries (detectors_tests PUBLIC detectors catch_common_${PROJECT_NAME} ${matio})
  add_dependencies (detectors_tests gen_data_q64_N60_w4 gen_data_q64_N60_w34 ${LIST_DETECTOR_DEP})

  catch_discover_tests (detectors_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # For generic detectors
  set (LIST_DETECTOR_TB "")
  set (LIST_DETECTOR_DEP "")
  foreach (TARGET_Q ${LIST_TARGET_Q})
    execute_process (COMMAND ${Python_EXECUTABLE} -c "import math as m; print(m.sqrt(${TARGET_Q}), end='')" OUTPUT_VARIABLE TARGET_SQRT_Q)

    set (LOCAL_FILENAME CDetectorSerialGeneric_q${TARGET_Q}_N60_w1_high_test.cpp)
    list (APPEND LIST_DETECTOR_TB ${LOCAL_FILENAME})
    configure_file (testing/CDetectorSerialGeneric_qXX_N60_w1_high_test.cpp.in ${LOCAL_FILENAME} @ONLY)
    unset (LOCAL_FILENAME)

    list (APPEND LIST_DETECTOR_DEP gen_data_q${TARGET_Q}_N60_w1)
    unset (TARGET_SQRT_Q)
  endforeach ()

  add_executable (
    detectors_generic_tests testing/CDetectorSerialGeneric_w1_low_test.cpp testing/CDetectorSerialGeneric_w4_pi_4_low_test.cpp
                            testing/CDetectorSerialGeneric_w34_pi_8_low_test.cpp ${LIST_DETECTOR_TB}
  )
  target_link_libraries (detectors_generic_tests PUBLIC detectors_generic catch_common_${PROJECT_NAME} ${matio})
  add_dependencies (detectors_generic_tests gen_data_q64_N60_w4 gen_data_q64_N60_w34 ${LIST_DETECTOR_DEP})

  catch_discover_tests (detectors_generic_tests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif ()
