#include <vector>

#include "./CCorrelationEngineGeneric.hpp"
#include "./TimeSliding/CCorrAbsMax/CCorrAbsMaxGeneric.hpp"
#include "./TimeSliding/CQSpannedSequentialAdder/CQSpannedSequentialAdderRC.hpp"

const std::vector<float> & QCSP::StandaloneDetector::CCorrelationEngineGeneric::get_pn() const { return corr_abs_max.get_pn(); }

float QCSP::StandaloneDetector::CCorrelationEngineGeneric::process(float re_in, float im_in) {
    float re_it, im_it;
    iterative_adder.process(re_in, im_in, &re_it, &im_it);
    return corr_abs_max.process(re_it, im_it);
}

QCSP::StandaloneDetector::CCorrelationEngineGeneric::CCorrelationEngineGeneric(const std::vector<float> & _pn)
    : iterative_adder(_pn.size()),
      corr_abs_max(_pn) {}
