#ifndef _QCSP_PASSED_DETECTION_STATE_TEMPLATE_HPP_
#define _QCSP_PASSED_DETECTION_STATE_TEMPLATE_HPP_ 1

#include <cstdint>

namespace QCSP {
namespace StandaloneDetector {

template <typename TScore, typename TFreq, unsigned Tp_omega>
struct DetectionState {
    static_assert(Tp_omega > 0U, "Can't have no hypothesis.");
    bool frame_detected;
    bool max_found;

    TScore scores[Tp_omega];
    TScore max_score;

    uint64_t chip_since_last_det;
    uint64_t chip_from_max;
    TFreq    frequency_offset; // The frequency offset corresponding to the maximum score
    uint32_t frequency_index;  // The index of the frequency offset in the frequency array
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_DETECTION_STATE_TEMPLATE_HPP_
