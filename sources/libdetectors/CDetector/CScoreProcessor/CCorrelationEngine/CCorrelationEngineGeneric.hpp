#ifndef _QCSP_PASSED_CORRELATION_GENERIC_HPP_
#define _QCSP_PASSED_CORRELATION_GENERIC_HPP_ 1

#include <vector>

#include "./TimeSliding/CCorrAbsMax/CCorrAbsMaxGeneric.hpp"
#include "./TimeSliding/CIterativeAdder/CIterativeAdderGeneric.hpp"

namespace QCSP {
namespace StandaloneDetector {

class CCorrelationEngineGeneric {
private:
    CIterativeAdderGeneric iterative_adder;
    CCorrAbsMaxGeneric     corr_abs_max;

public:
    const std::vector<float> & get_pn() const { return corr_abs_max.get_pn(); }

    virtual float process(float re_in, float im_in);

    CCorrelationEngineGeneric(const std::vector<float> & pn);

    virtual ~CCorrelationEngineGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORRELATION_GENERIC_HPP_
