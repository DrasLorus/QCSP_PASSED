#include "./CDetector.hpp"

using namespace QCSP::StandaloneDetector;

namespace {
class phantom : public CDetectorSerial<60, 64, 9, float, true> {};
} // namespace
