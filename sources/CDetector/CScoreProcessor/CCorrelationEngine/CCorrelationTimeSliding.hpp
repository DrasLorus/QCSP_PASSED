#ifndef _QCSP_PASSED_CORRELATION_TIME_SLIDING_HPP_
#define _QCSP_PASSED_CORRELATION_TIME_SLIDING_HPP_ 1

#include "./TimeSliding/CCorrAbsMax/CCorrAbsMax.hpp"
#include "./TimeSliding/CIterativeAdder/CIterativeAdder.hpp"
#include "CCorrelationEngineTemplate.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned Tq, typename TIn_Type = float>
class CCorrelationTimeSliding : public CCorrelationEngineTemplate<Tq, TIn_Type> {
    using CCorrelationEngineTemplate<Tq, TIn_Type>::q;

private:
    CIterativeAdder<q, TIn_Type> iterative_adder;
    CCorrAbsMax<q, TIn_Type>     corr_abs_max;

public:
    virtual const TIn_Type * get_pn() const override { return corr_abs_max.get_pn(); }

    virtual TIn_Type process(TIn_Type re_in, TIn_Type im_in) override {
        TIn_Type re_it, im_it;
        iterative_adder.process(re_in, im_in, &re_it, &im_it);
        return corr_abs_max.process(re_it, im_it);
    }

    template <typename Tpn>
    CCorrelationTimeSliding(Tpn * pn)
        : CCorrelationEngineTemplate<Tq, TIn_Type>(),
          iterative_adder(),
          corr_abs_max(pn) {}

    virtual ~CCorrelationTimeSliding() = default;
};

template <unsigned Tq>
class CCorrelationTimeSliding<Tq, int16_t> : public CCorrelationEngineTemplate<Tq, int16_t> {
    using CCorrelationEngineTemplate<Tq, int16_t>::q;

private:
    CIterativeAdder<q, int16_t> iterative_adder;
    CCorrAbsMax<q, int16_t>     corr_abs_max;

public:
#ifndef USE_PN_XOR_TRICK
    virtual const int8_t * get_pn() const override { return corr_abs_max.get_pn(); }
#else
    virtual const bool * get_pn() const override { return corr_abs_max.get_pn(); }
#endif

    virtual uint32_t process(int16_t re_in, int16_t im_in) override {
        int32_t re_it, im_it;
        iterative_adder.process(re_in, im_in, &re_it, &im_it);
        return corr_abs_max.process(re_it, im_it);
    }

    template <typename Tpn>
    CCorrelationTimeSliding(Tpn * pn)
        : CCorrelationEngineTemplate<Tq, int16_t>(),
          iterative_adder(),
          corr_abs_max(pn) {}

    virtual ~CCorrelationTimeSliding() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORRELATION_TIME_SLIDING_HPP_
