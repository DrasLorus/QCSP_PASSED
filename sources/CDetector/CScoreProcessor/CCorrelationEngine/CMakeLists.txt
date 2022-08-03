add_subdirectory (TimeSliding/CIterativeAdder)

add_subdirectory (TimeSliding/CCorrAbsMax)

add_library (correngine OBJECT CCorrelationEngine.cpp)

target_link_libraries (correngine PUBLIC iterativeadder corrabsmax)

if (ENABLE_TESTING)
  include (CTest)
  include (Catch)

  # TODO Add tests

endif ()
