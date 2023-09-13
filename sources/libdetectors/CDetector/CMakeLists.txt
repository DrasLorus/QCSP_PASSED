add_subdirectory (CScoreProcessor)

add_library (detectors STATIC CDetector.cpp)
target_include_directories (detectors PUBLIC ${LIB_ROOT_DIRECTORY})
target_link_libraries (detectors PUBLIC score_processor)

add_library (detectors_generic STATIC CDetectorSerialGeneric.cpp CDetectionStateGeneric.cpp)
target_include_directories (detectors_generic PUBLIC ${LIB_ROOT_DIRECTORY})
target_link_libraries (detectors_generic PUBLIC score_processor_generic)

list (
  APPEND
  DETECTORS_GENERIC_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetectionStateGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetectorSerialGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetectorGenericInterface.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CCorrelationEngine/CCorrelationEngineGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CCorrelationEngine/TimeSliding/CCorrAbsMax/CCorrAbsMaxGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CCorrelationEngine/TimeSliding/CIterativeAdder/CIterativeAdderGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CNorm/CNormGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CScoreAccumulator/CScoreAccumulatorGeneric.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CScoreProcessor/CScoreProcessorGeneric.hpp
)

list (
  APPEND
  DETECTORS_GENERIC_INCLUDE_HEADERS
  ${CMAKE_SOURCE_DIR}/includes/CDetectionStateGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CDetectorSerialGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CDetectorGenericInterface.hpp
  ${CMAKE_SOURCE_DIR}/includes/CCorrelationEngineGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CCorrAbsMaxGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CIterativeAdderGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CNormGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CScoreAccumulatorGeneric.hpp
  ${CMAKE_SOURCE_DIR}/includes/CScoreProcessorGeneric.hpp
)

file (MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/includes)
foreach (DG_HEADERS IN ZIP_LISTS DETECTORS_GENERIC_HEADERS DETECTORS_GENERIC_INCLUDE_HEADERS)
  if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
    find_program (POWERSHELL_PATH NAMES powershell REQUIRED)
    add_custom_command (
      OUTPUT ${DG_HEADERS_1} COMMAND ${POWERSHELL_PATH} -Command
                                     "(cat ${DG_HEADERS_0})  -creplace '#include \"\\./.*/','#include \"' > ${DG_HEADERS_1}"
    )
  else ()
    find_program (SED_PATH NAMES sed REQUIRED)
    add_custom_command (
      OUTPUT ${DG_HEADERS_1} COMMAND sh ${CMAKE_SOURCE_DIR}/scripts/path-remover.sh ${SED_PATH} ${DG_HEADERS_0} ${DG_HEADERS_1}
    ) # REVIEW
  endif ()
endforeach ()
add_custom_target (detectors_generic_public_headers DEPENDS ${DETECTORS_GENERIC_INCLUDE_HEADERS})

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
