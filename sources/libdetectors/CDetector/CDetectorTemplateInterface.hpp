#ifndef _QCSP_PASSED_DETECTOR_TEMPLATE_INTERFACE_HPP_
#define _QCSP_PASSED_DETECTOR_TEMPLATE_INTERFACE_HPP_ 1

#include <stdexcept>
#include <string>
#include <vector>

#include "./CDetectionStateTemplate.hpp"
#include "./CScoreProcessor/CScoreProcessor.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned TFrameSize, unsigned Tq, unsigned Tp_omega, typename TIn_Type = float, bool normed = true, CorrelationEngineType variant = TIME_SLIDING>
class CDetector {
public:
    using state_t = DetectionState<TIn_Type, TIn_Type, Tp_omega>;

    virtual void process(TIn_Type re_in, TIn_Type im_in, state_t * state)     = 0;
    virtual void process_sqr(TIn_Type re_in, TIn_Type im_in, state_t * state) = 0;

    virtual ~CDetector() = default;
};

// 16-bits integer specialization
template <unsigned TFrameSize, unsigned Tq, unsigned Tp_omega, bool normed, CorrelationEngineType variant>
class CDetector<TFrameSize, Tq, Tp_omega, int16_t, normed, variant> {
public:
    using state_t = DetectionState<uint32_t, float, Tp_omega>;

    virtual void process(int16_t re_in, int16_t im_in, state_t * state)     = 0;
    virtual void process_sqr(int16_t re_in, int16_t im_in, state_t * state) = 0;

    virtual ~CDetector() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif
