#ifndef _QCSP_PASSED_NORM_GENERIC_HPP_
#define _QCSP_PASSED_NORM_GENERIC_HPP_ 1

#include <cstdint>
#include <vector>

namespace QCSP {
namespace StandaloneDetector {

/**
 * @brief Compute the square value of 2-Norm of a q-long complex vector
 *
 * @details The norm is computed iteratively on magnitudes, so the square value comes naturally.
 *
 * @tparam Tq size of the vector
 * @tparam float data representation
 */
class CNormGeneric {
public:
    const unsigned q;    // = Tq;
    const unsigned mask; // = Tq - 1;

private:
    std::vector<float> magn_fifo; // [q]
    float              norm_accumulator;

    uint32_t counter;

public:
    virtual float process(float re_in, float im_in);

    CNormGeneric(unsigned q);

    virtual ~CNormGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_NORM_HPP_
