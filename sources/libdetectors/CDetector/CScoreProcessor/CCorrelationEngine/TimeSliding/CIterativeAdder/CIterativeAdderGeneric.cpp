#include <cstring>
#include <stdexcept>
#include <vector>

#include "./CIterativeAdderGeneric.hpp"

void QCSP::StandaloneDetector::CIterativeAdderGeneric::process(float re_in, float im_in, float * re_out, float * im_out) {
    const unsigned curr_cnt = counter;

    const float old_re = fifo_re[curr_cnt];
    const float old_im = fifo_im[curr_cnt];

    fifo_re[curr_cnt] = re_in;
    fifo_im[curr_cnt] = im_in;

    counter = (curr_cnt + 1) & cnt_mask;

    *re_out = re_in - old_re;
    *im_out = im_in - old_im;
}

QCSP::StandaloneDetector::CIterativeAdderGeneric::CIterativeAdderGeneric(unsigned _q)
    : q(_q),
      cnt_mask(_q - 1),
      fifo_re(_q, 0),
      fifo_im(_q, 0) {
    if (_q & 1) {
        throw std::runtime_error("q must be a power of 2.");
    }
    counter = 0;
}
