#ifndef _QCSP_PASSED_NORM_HPP_
#define _QCSP_PASSED_NORM_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/metatypes.hpp"
#include "Miscellanous/misc.hpp"
#include <type_traits>

namespace QCSP {
namespace StandaloneDetector {

/**
 * @brief Compute the square value of 2-Norm of a q-long complex vector
 *
 * @details The norm is computed iteratively on magnitudes, so the square value comes naturally.
 *
 * @tparam Tq size of the vector
 * @tparam Tin_type data representation
 */
template <unsigned Tq, typename Tin_type = float,
          int Tinput_width = sizeof(Tin_type) * 8, int Tinput_int = sizeof(Tin_type) * 4,
          bool Tis_float = std::is_floating_point<Tin_type>::value>
class CNorm {
    static_assert(is_pow2(Tq), "q must be a power of 2.");
    static_assert(std::is_floating_point<Tin_type>::value, "template parameters overspecified or Tin_type is neither a floating-point nor an integral type");

public:
    static constexpr unsigned q    = Tq;
    static constexpr unsigned mask = Tq - 1;

private:
    Tin_type magn_fifo[q];
    Tin_type norm_accumulator;

    uint32_t counter;

public:
    Tin_type process(Tin_type re_in, Tin_type im_in) {
        const uint32_t curr_counter = counter;

        const Tin_type old_magn  = magn_fifo[curr_counter];
        const Tin_type curr_sqrn = norm_accumulator;

        const Tin_type new_magn = re_in * re_in + im_in * im_in;

        const Tin_type new_sqrn = curr_sqrn + new_magn - old_magn;

        norm_accumulator        = new_sqrn;
        counter                 = (counter + 1) & mask;
        magn_fifo[curr_counter] = new_magn;

        return new_sqrn;
    }

    CNorm()
        : norm_accumulator(Tin_type(1)),
          counter(0) {
        memset(magn_fifo, 0, sizeof(Tin_type) * q);
        magn_fifo[q - 1] = Tin_type(1);
    }

    virtual ~CNorm() = default;
};

/**
 * @brief CNorm Specialization for 16 bit input.
 *
 * @details Compute the square value of 2-Norm of a q-long complex vector.
 *  The norm is computed iteratively on magnitudes, so the square value comes naturally.
 *  This specialization assume a 16-bit quantified input, as delivered by USRPs
 *
 * @tparam Tq size of the vector
 */
template <unsigned Tq, typename Tin_type, int Tinput_width, int Tinput_int>
class CNorm<Tq, Tin_type, Tinput_width, Tinput_int, false> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");
    static_assert(std::is_integral<Tin_type>::value, "template parameters overspecified or Tin_type is neither a floating-point nor an integral type");

public:
    static constexpr unsigned q          = Tq;
    static constexpr unsigned p          = pow2_log2<q>();
    static constexpr unsigned mask       = Tq - 1;
    static constexpr unsigned in_width   = Tinput_width;
    static constexpr int64_t  in_int     = Tinput_int;
    static constexpr int64_t  in_quote   = Tinput_width - Tinput_int;
    static constexpr unsigned fifo_width = 2 * in_width + 1 - in_quote;
    static constexpr unsigned fifo_int   = 2 * in_int + 1;
    static constexpr unsigned reg_width  = 2 * in_width + 2 * (p + 1) - in_quote;
    static constexpr unsigned reg_int    = 2 * in_int + 2 * (p + 1);

    using input_t      = intxx_t<in_width, metatypes::FROM_BITS>;
    using fifo_t       = uintxx_t<fifo_width, metatypes::FROM_BITS>;
    using local_accu_t = uintxx_t<reg_width + 2, metatypes::FROM_BITS>;
    using accu_t       = uintxx_t<reg_width, metatypes::FROM_BITS>;
    using counter_t    = uintxx_t<q, metatypes::FROM_MAX>;
    using prod_t       = uintxx_t<2 * in_width, metatypes::FROM_BITS>;
    using magn_full_t  = uintxx_t<2 * in_width + 1, metatypes::FROM_BITS>;
    using magn_t       = uintxx_t<2 * in_width + 1 - in_quote, metatypes::FROM_BITS>;
    using out_t        = typename truncation<reg_width, p + 1>::out_t;

    static constexpr magn_t initial_value = 1ULL << (in_int * 2 + 1);
    static constexpr out_t  initial_out   = truncation<reg_width, p + 1>::apply(initial_value);

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

        const magn_t new_magn = magn_t(magn_full_t(prod_t(re_in * re_in) + prod_t(im_in * im_in)) >> in_quote);

        const local_accu_t full_accu      = local_accu_t(curr_sqrn + new_magn) - old_magn;
        const local_accu_t truncated_accu = truncation<reg_width + 2, 1, metatypes::UNSIGNED, false>::apply(full_accu);
        const accu_t       new_sqrn       = saturation<reg_width + 1, 1>::apply(truncated_accu);

        norm_accumulator        = new_sqrn;
        magn_fifo[curr_counter] = new_magn;

        const out_t trunc_sqrn = truncation<reg_width, p + 1>::apply(new_sqrn);

        //* To deeply debug, include <cstdio> and uncomment the following
        // fprintf(stderr, "%-8u %-8u %-8lu %-8u %-8lu %-8lu %-8u\n", re_in, im_in, old_magn, curr_counter, new_magn, new_sqrn, trunc_sqrn);

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

template <unsigned Tq>
using CNorm_i16 = CNorm<Tq, int16_t, 16, 4>;

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_NORM_HPP_
