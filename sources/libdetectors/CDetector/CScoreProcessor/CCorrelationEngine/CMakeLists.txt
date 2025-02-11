add_subdirectory (TimeSliding/CQSpannedSequentialAdder)

add_subdirectory (TimeSliding/CCorrAbsMax)

add_library (correngine STATIC CCorrelationEngine.cpp)
target_link_libraries (correngine PUBLIC qspannedsequentialadder corrabsmax)

add_library (correngine_generic STATIC CCorrelationEngineGeneric.cpp)
target_link_libraries (correngine_generic PUBLIC qspannedsequentialadder_rc corrabsmax_generic)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  # TODO Add tests

endif ()
