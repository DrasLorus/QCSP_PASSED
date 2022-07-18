#include "CNorm.hpp"

namespace {
class phantom : public QCSP::StandaloneDetector::CNorm<64, float> {};
} // namespace

namespace {
class phantom_int : public QCSP::StandaloneDetector::CNorm<64, int> {};
} // namespace
