#ifndef _QCSP_PASSED_DETECTOR_GENERIC_INTERFACE_HPP_
#define _QCSP_PASSED_DETECTOR_GENERIC_INTERFACE_HPP_ 1

#include <vector>

#include "./CDetectionStateGeneric.hpp"

namespace QCSP {
namespace StandaloneDetector {

class CDetectorGeneric {
public:
    using state_t = DetectionStateGeneric;

    virtual void process(float re_in, float im_in, state_t * state)     = 0;
    virtual void process_sqr(float re_in, float im_in, state_t * state) = 0;

    virtual ~CDetectorGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_DETECTOR_GENERIC_INTERFACE_HPP_
