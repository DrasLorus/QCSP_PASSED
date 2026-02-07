#ifndef _QCSP_PASSED_DETECTION_STATE_GENERIC_HPP_
#define _QCSP_PASSED_DETECTION_STATE_GENERIC_HPP_ 1

#include <cstdint>
#include <vector>

namespace QCSP {
namespace StandaloneDetector {

struct DetectionStateGeneric {
    const uint16_t p_omega;

    bool frame_detected;
    bool max_found;

    std::vector<float> scores;
    float              max_score;

    uint64_t chip_since_last_det;
    uint64_t chip_from_max;
    float    frequency_offset; // The frequency offset corresponding to the maximum score
    uint32_t frequency_index;  // The index of the frequency offset in the frequency array

    DetectionStateGeneric() = delete;
    DetectionStateGeneric(uint16_t _p_omega);
    ~DetectionStateGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_DETECTION_STATE_GENERIC_HPP_
