#include "CNorm.hpp"
#include "CNormFP.hpp"

namespace {
class phantom : public QCSP::StandaloneDetector::CNorm<64, float> {};
} // namespace

namespace {
class phantom_int : public QCSP::StandaloneDetector::CNorm<64, int16_t> {};
} // namespace
