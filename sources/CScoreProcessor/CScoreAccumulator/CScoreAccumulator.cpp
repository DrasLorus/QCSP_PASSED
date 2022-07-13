#include "CScoreAccumulator.hpp"

namespace {
class phantom : public QCSP::StandaloneDetector::CScoreAccumulator<60, 64, float> {};
} // namespace

namespace {
class phantom_int : public QCSP::StandaloneDetector::CScoreAccumulator<60, 64, uint32_t> {};
} // namespace
