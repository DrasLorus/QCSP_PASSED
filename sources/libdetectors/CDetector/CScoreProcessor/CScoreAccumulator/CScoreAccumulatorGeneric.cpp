#include <vector>

#include "./CScoreAccumulatorGeneric.hpp"

float QCSP::StandaloneDetector::CScoreAccumulatorGeneric::process(float new_max) {
    const uint32_t curr_counter  = counter;
    const uint32_t score_counter = curr_counter & mask;
    const float    old_max       = fifos_max[curr_counter];
    const float    old_score     = score_registers[score_counter];

    counter = (curr_counter + 1) * uint32_t(curr_counter != (N * q - 1)); // Faster !
    // counter = (curr_counter + 1) % (N * q);

    const float new_score = old_score + new_max - old_max;

    score_registers[score_counter] = new_score;
    fifos_max[curr_counter]        = new_max;

    return new_score;
}

QCSP::StandaloneDetector::CScoreAccumulatorGeneric::CScoreAccumulatorGeneric(unsigned _q, unsigned _N)
    : N(_N),
      q(_q),
      mask(q - 1),
      fifos_max(N * q, 0),
      score_registers(q, 0),
      counter(0) {
}
