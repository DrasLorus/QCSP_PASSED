#include "./CDetectionStateGeneric.hpp"

QCSP::StandaloneDetector::DetectionStateGeneric::DetectionStateGeneric(uint16_t _p_omega)
    : p_omega(_p_omega),
      frame_detected(false),
      max_found(false),
      scores(_p_omega, 0),
      max_score(0.0f),
      chip_since_last_det(0),
      chip_from_max(0),
      frequency_offset(0.0f), // The frequency offset corresponding to the maximum score
      frequency_index(0) {}
