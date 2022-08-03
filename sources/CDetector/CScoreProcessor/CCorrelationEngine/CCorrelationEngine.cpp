#include "./CCorrelationEngine.hpp"
#include "./CCorrelationTimeSliding.hpp"

namespace {
class phantom : public QCSP::StandaloneDetector::CCorrelationTimeSliding<64, float> {};
} // namespace