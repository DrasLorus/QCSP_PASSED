#include <cstdint>
#include <vector>

#include "./CNormGeneric.hpp"

float QCSP::StandaloneDetector::CNormGeneric::process(float re_in, float im_in) {
    const uint32_t curr_counter = counter;

    const float old_magn  = magn_fifo[curr_counter];
    const float curr_sqrn = norm_accumulator;

    const float new_magn = re_in * re_in + im_in * im_in;

    const float new_sqrn = curr_sqrn + new_magn - old_magn;

    norm_accumulator        = new_sqrn;
    counter                 = (counter + 1) & mask;
    magn_fifo[curr_counter] = new_magn;

    return new_sqrn;
}

QCSP::StandaloneDetector::CNormGeneric::CNormGeneric(unsigned _q)
    : q(_q),
      mask(q - 1),
      magn_fifo(q, 0),
      norm_accumulator(float(1)),
      counter(0) {
    magn_fifo[q - 1] = float(1);
}
