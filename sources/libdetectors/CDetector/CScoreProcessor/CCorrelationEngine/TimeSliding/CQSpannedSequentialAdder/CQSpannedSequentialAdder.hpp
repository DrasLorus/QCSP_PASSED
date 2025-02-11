#ifndef _QCSP_PASSED_ITERATIVE_ADDER_HPP_
#define _QCSP_PASSED_ITERATIVE_ADDER_HPP_

#include <cstring>

#include "Miscellanous/misc.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned q, typename Tin_type = float>
class CQSpannedSequentialAdder {
    static_assert(is_pow2(q), "q must be a power of 2!");

public:
    static constexpr unsigned cnt_mask = q - 1;

private:
    Tin_type fifo_re[q];
    Tin_type fifo_im[q];
    unsigned counter;

public:
    void process(Tin_type re_in, Tin_type im_in, Tin_type * re_out, Tin_type * im_out) {
        const unsigned curr_cnt = counter;

        const Tin_type old_re = fifo_re[curr_cnt];
        const Tin_type old_im = fifo_im[curr_cnt];

        fifo_re[curr_cnt] = re_in;
        fifo_im[curr_cnt] = im_in;

        counter = (curr_cnt + 1) & cnt_mask;

        *re_out = re_in - old_re;
        *im_out = im_in - old_im;
    }

    CQSpannedSequentialAdder() {
        memset(fifo_re, 0, sizeof(Tin_type) * q);
        memset(fifo_im, 0, sizeof(Tin_type) * q);
        counter = 0;
    }

    ~CQSpannedSequentialAdder() = default;
};

/**
 * @brief CQSpannedSequentialAdder Specialization for 16 bit input
 *
 * @details This specialization assume a 16-bit quantified input, as delivered by USRPs
 *
 * @tparam q size of the QCSP sequence
 */
template <unsigned q>
class CQSpannedSequentialAdder<q, int16_t> {
    static_assert(is_pow2(q), "q must be a power of 2!");

public:
    static constexpr unsigned cnt_mask = q - 1;

private:
    int16_t  fifo_re[q];
    int16_t  fifo_im[q];
    unsigned counter;

public:
    void process(int16_t re_in, int16_t im_in, int32_t * re_out, int32_t * im_out) {
        const unsigned curr_cnt = counter;

        const int16_t old_re = fifo_re[counter];
        const int16_t old_im = fifo_im[counter];

        fifo_re[counter] = re_in;
        fifo_im[counter] = im_in;

        counter = (curr_cnt + 1) & cnt_mask;

        *re_out = re_in - old_re;
        *im_out = im_in - old_im;
    }

    CQSpannedSequentialAdder() : counter(0) {
        memset(fifo_re, 0, sizeof(int16_t) * q);
        memset(fifo_im, 0, sizeof(int16_t) * q);
    }

    ~CQSpannedSequentialAdder() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_ITERATIVE_ADDER_HPP_
