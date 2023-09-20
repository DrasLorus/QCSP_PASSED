set (LIB_ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (CDetector)

list (
  APPEND
  DETECTORS_GENERIC_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CDetectionStateGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CDetectorSerialGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CCorrelationEngine/CCorrelationEngineGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CCorrelationEngine/TimeSliding/CCorrAbsMax/CCorrAbsMaxGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CCorrelationEngine/TimeSliding/CIterativeAdder/CIterativeAdderGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CNorm/CNormGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CScoreAccumulator/CScoreAccumulatorGeneric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CDetector/CScoreProcessor/CScoreProcessorGeneric.cpp
)

add_library (qcsp_passed STATIC ${DETECTORS_GENERIC_SOURCES})
set_target_properties (qcsp_passed PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
