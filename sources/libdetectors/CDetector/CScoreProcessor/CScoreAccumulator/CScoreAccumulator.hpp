#include <cstdio>
#ifndef _QCSP_PASSED_SCORE_ACCUMULATOR_HPP_
#define _QCSP_PASSED_SCORE_ACCUMULATOR_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned TFrameSize, unsigned Tq, class Tin_type = float>
class CScoreAccumulator {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

private:
    static constexpr uint32_t N    = TFrameSize;
    static constexpr uint32_t q    = Tq;
    static constexpr uint32_t mask = q - 1;

    Tin_type fifos_max[N * q];   // linearized max fifos
    Tin_type score_registers[q]; // aggregated score registers

    uint32_t counter;

public:
    Tin_type process(Tin_type new_max) {
        const uint32_t curr_counter  = counter;
        const uint32_t score_counter = curr_counter & mask;
        const Tin_type old_max       = fifos_max[curr_counter];
        const Tin_type old_score     = score_registers[score_counter];

        counter = (curr_counter + 1) * uint32_t(curr_counter != (N * q - 1)); // Faster !
        // counter = (curr_counter + 1) % (N * q);

        const Tin_type new_score = old_score + new_max - old_max;

        score_registers[score_counter] = new_score;
        fifos_max[curr_counter]        = new_max;

        return new_score;
    }

    CScoreAccumulator() : counter(0) {
        memset(fifos_max, 0, N * q * sizeof(Tin_type));
        memset(score_registers, 0, q * sizeof(Tin_type));
    }

    virtual ~CScoreAccumulator() = default;
};

/**
 * @brief CScoreAccumulator Specialization for 16 bit input in the score processor
 *
 * @details This specialization assume a 16-bit quantified input for the full score processor, as delivered by USRPs
 *          It assumes the corr_abs_max value is quantified on 24 bits.
 *
 * @tparam TFrameSize size of the QCSP frame
 * @tparam Tq size of the QCSP sequence
 */
template <unsigned TFrameSize, unsigned Tq>
class CScoreAccumulator<TFrameSize, Tq, int16_t> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr uint32_t N    = TFrameSize;
    static constexpr uint32_t q    = Tq;
    static constexpr uint32_t mask = q - 1;

private:
    uint32_t fifos_max[N * q];   // linearized max fifos (24 bits)
    uint64_t score_registers[q]; // aggegated score registers (30 bits)

    uint32_t counter;

public:
    uint32_t process(uint32_t new_max) {
        const uint32_t curr_counter  = counter;
        const uint32_t score_counter = curr_counter & mask;
        const uint32_t old_max       = fifos_max[curr_counter];
        const uint32_t old_score     = score_registers[score_counter];

        counter = (curr_counter + 1) * uint32_t(curr_counter != (N * q - 1)); // Faster !
        // counter = (curr_counter + 1) % (N * q);

        const uint32_t new_score = (old_score + new_max) - old_max;

        score_registers[score_counter] = new_score;
        fifos_max[curr_counter]        = new_max;

        return new_score;
    }

    CScoreAccumulator() : counter(0) {
        memset(fifos_max, 0, N * q * sizeof(uint32_t));
        memset(score_registers, 0, q * sizeof(uint64_t));
    }

    virtual ~CScoreAccumulator() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_SCORE_ACCUMULATOR_HPP_
