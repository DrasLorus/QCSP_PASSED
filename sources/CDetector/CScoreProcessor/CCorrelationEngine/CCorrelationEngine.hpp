#ifndef _QCSP_PASSED_CORRELATION_ENGINE_HPP_
#define _QCSP_PASSED_CORRELATION_ENGINE_HPP_ 1

#include "./CCorrelationEngineTemplate.hpp"
#include "./CCorrelationTimeSliding.hpp"

namespace QCSP {
namespace StandaloneDetector {

enum CorrelationEngineType : int {
    TIME_SLIDING = 0,
    FFT          = 1
};

template <unsigned Tq, typename TIn_Type = float, CorrelationEngineType TType = TIME_SLIDING>
class CCorrelationEngine;

template <unsigned Tq, typename TIn_Type>
class CCorrelationEngine<Tq, TIn_Type, TIME_SLIDING> : public CCorrelationTimeSliding<Tq, TIn_Type> {
public:
    template <typename Tpn>
    CCorrelationEngine(Tpn * pn) : CCorrelationTimeSliding<Tq, TIn_Type>(pn) {}
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORRELATION_ENGINE_HPP_