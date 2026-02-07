#ifndef _QCSP_PASSED_ITERATIVE_ADDER_GENERIC_HPP_
#define _QCSP_PASSED_ITERATIVE_ADDER_GENERIC_HPP_ 1

#include <vector>

namespace QCSP {
namespace StandaloneDetector {

class CQSpannedSequentialAdderRC {
public:
    const unsigned q;
    const unsigned cnt_mask; // = q - 1;

private:
    std::vector<float> fifo_re; // [q]
    std::vector<float> fifo_im; // [q]
    unsigned           counter;

public:
    virtual void process(float re_in, float im_in, float * re_out, float * im_out);

    CQSpannedSequentialAdderRC(unsigned q);

    virtual ~CQSpannedSequentialAdderRC() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_ITERATIVE_ADDER_GENERIC_HPP_
