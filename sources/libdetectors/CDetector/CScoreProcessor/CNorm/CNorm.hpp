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
    static_assert(std::is_floating_point<Tin_type>::value, "wrong template parameters overspecified or Tin_type is not a floating-point type");

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
} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_NORM_HPP_
