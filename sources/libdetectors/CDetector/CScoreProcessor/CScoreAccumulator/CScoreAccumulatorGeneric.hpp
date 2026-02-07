#ifndef _QCSP_PASSED_SCORE_ACCUMULATOR_GENERIC_HPP_
#define _QCSP_PASSED_SCORE_ACCUMULATOR_GENERIC_HPP_ 1

#include <cstdint>
#include <vector>

namespace QCSP {
namespace StandaloneDetector {

class CScoreAccumulatorGeneric {

private:
    const uint32_t N;    // = TFrameSize;
    const uint32_t q;    // = Tq;
    const uint32_t mask; // = q - 1;

    std::vector<float> fifos_max;       // [N * q] - linearized max fifos
    std::vector<float> score_registers; // [q]     - aggregated score registers

    uint32_t counter;

public:
    virtual float process(float new_max);

    CScoreAccumulatorGeneric(unsigned q, unsigned N);

    virtual ~CScoreAccumulatorGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_SCORE_ACCUMULATOR_GENERIC_HPP_
