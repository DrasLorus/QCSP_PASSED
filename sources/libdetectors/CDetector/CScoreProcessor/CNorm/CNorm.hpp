#ifndef _QCSP_PASSED_NORM_HPP_
#define _QCSP_PASSED_NORM_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

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
template <unsigned Tq, typename Tin_type = float>
class CNorm {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

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
template <unsigned Tq>
class CNorm<Tq, int16_t> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q    = Tq;
    static constexpr unsigned mask = Tq - 1;

private:
    uint64_t magn_fifo[q];
    uint64_t norm_accumulator;

    uint32_t counter;

public:
    uint16_t process_sqr(int16_t re_in, int16_t im_in) {
        const uint32_t curr_counter = counter;

        const uint64_t curr_sqrn = norm_accumulator;
        const uint64_t old_magn  = magn_fifo[curr_counter];

        // 2 * 16 + 1 > 32 so uint64_t
        const uint64_t new_magn = uint32_t(re_in * re_in) + uint32_t(im_in * im_in);

        const uint64_t new_sqrn = (curr_sqrn + new_magn) - old_magn;

        norm_accumulator        = new_sqrn;
        counter                 = (curr_counter + 1) & mask;
        magn_fifo[curr_counter] = new_magn;

        const uint16_t trunc_sqrn = uint16_t(new_sqrn >> (16 + pow2_log2<q>() + 1));

        //* To deeply debug, include <cstdio> and uncomment the following
        // fprintf(stderr, "%-8u %-8u %-8lu %-8u %-8lu %-8lu %-8u\n", re_in, im_in, old_magn, curr_counter, new_magn, new_sqrn, trunc_sqrn);

        return trunc_sqrn;
    }

    CNorm()
        : norm_accumulator(1LLU << (16 + pow2_log2<q>() + 1)),
          counter(0) {
        memset(magn_fifo, 0, q * sizeof(uint64_t));
        magn_fifo[q - 1] = 1LLU << (16 + pow2_log2<q>() + 1);
    }

    virtual ~CNorm() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_NORM_HPP_
