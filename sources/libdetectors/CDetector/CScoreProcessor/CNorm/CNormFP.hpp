#ifndef _QCSP_PASSED_NORM_FP_HPP_
#define _QCSP_PASSED_NORM_FP_HPP_ 1

#include "CNorm.hpp"

#include "Miscellanous/metatypes.hpp"
#include "Miscellanous/misc.hpp"
#include <type_traits>

namespace QCSP {
namespace StandaloneDetector {

/**
 * @brief CNorm Specialization for fixed-point or integer input.
 *
 * @details Compute the square value of 2-Norm of a q-long complex vector.
 *  The norm is computed iteratively on magnitudes, so the square value comes naturally.
 *  This specialization does not assume any quantization of the input.
 *  Quantization are deduced during compilation.
 *
 * @tparam Tq size of the vector
 * @tparam Tinput_width size of the input (in bits)
 * @tparam Tinput_int size of the input dedicated to the integer part (in bits)
 */
template <unsigned Tq, int Tinput_width, int Tinput_int>
class CNorm<Tq, intxx_t<Tinput_width, metatypes::FROM_BITS>, Tinput_width, Tinput_int, false> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q          = Tq;
    static constexpr unsigned p          = pow2_log2<q>();
    static constexpr unsigned mask       = Tq - 1;
    static constexpr unsigned in_width   = Tinput_width;
    static constexpr int64_t  in_int     = Tinput_int;
    static constexpr int64_t  in_quote   = Tinput_width - Tinput_int;
    static constexpr unsigned fifo_width = in_width;
    static constexpr unsigned fifo_int   = 2 * in_int;
    static constexpr unsigned reg_width  = fifo_width + p + 1;
    static constexpr unsigned reg_int    = 2 * in_int + p + 1;

    using out_stage = saturation<reg_width, p + 1>;

    using input_t      = intxx_t<in_width, metatypes::FROM_BITS>;
    using fifo_t       = uintxx_t<fifo_width, metatypes::FROM_BITS>;
    using local_accu_t = uintxx_t<reg_width + 2, metatypes::FROM_BITS>;
    using accu_t       = uintxx_t<reg_width, metatypes::FROM_BITS>;
    using counter_t    = uintxx_t<q, metatypes::FROM_MAX>;
    using prod_t       = uintxx_t<2 * in_width, metatypes::FROM_BITS>;
    using magn_full_t  = uintxx_t<2 * in_width + 1, metatypes::FROM_BITS>;
    using magn_t       = uintxx_t<2 * in_width + 1 - in_quote, metatypes::FROM_BITS>;
    using out_t        = typename out_stage::out_t;

    static constexpr magn_t initial_value = 1ULL << (in_int * 2 + 1);
    static constexpr out_t  initial_out   = out_stage::apply(initial_value);

private:
    fifo_t magn_fifo[q];
    accu_t norm_accumulator;

    counter_t counter;

public:
    out_t process_sqr(input_t re_in, input_t im_in) {
        const counter_t curr_counter = counter;
        counter                      = (curr_counter + 1) & mask;

        const accu_t curr_sqrn = norm_accumulator;
        const fifo_t old_magn  = magn_fifo[curr_counter];

        const auto   full_magn = magn_full_t(prod_t(re_in * re_in) + prod_t(im_in * im_in));
        const magn_t new_magn  = truncation<2 * in_width + 1, abs(in_width), metatypes::UNSIGNED>::apply(full_magn);

        const local_accu_t full_accu      = local_accu_t(curr_sqrn + new_magn) - old_magn;
        const local_accu_t truncated_accu = truncation<reg_width + 2, 1, metatypes::UNSIGNED, false>::apply(full_accu);
        const accu_t       new_sqrn       = saturation<reg_width + 1, 1>::apply(truncated_accu);

        norm_accumulator        = new_sqrn;
        magn_fifo[curr_counter] = new_magn;

        const auto  trunc_sqrn_tmp = out_stage::apply(new_sqrn);
        const out_t trunc_sqrn     = trunc_sqrn_tmp;

        //* To deeply debug, include <cstdio> and uncomment the following
        // fprintf(stderr, "%-8d\t%-8d\t%-8u\t%-8u\t%-8u\t%-8u\t%-8u\n", re_in, im_in, old_magn, curr_counter, new_magn, new_sqrn, trunc_sqrn);

        return trunc_sqrn;
    }

    CNorm()
        : norm_accumulator(initial_value),
          counter(0) {
        memset(magn_fifo, 0, q * sizeof(fifo_t));
        magn_fifo[q - 1] = initial_value;
    }

    virtual ~CNorm() = default;
};

template <unsigned Tq, int Tinput_width, int Tinput_int>
CNorm<Tq, intxx_t<Tinput_width, metatypes::FROM_BITS>, Tinput_width, Tinput_int, false> * make_new_cnorm_fp() {
    return new CNorm<Tq, intxx_t<Tinput_width, metatypes::FROM_BITS>, Tinput_width, Tinput_int, false>();
}

template <unsigned Tq>
using CNorm_i16_m10dB = CNorm<Tq, int16_t, 16, 4>;

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_NORM_FP_HPP_
