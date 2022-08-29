#ifndef _QCSP_PASSED_CORRELATION_ENGINE_TEMPLATE_HPP_
#define _QCSP_PASSED_CORRELATION_ENGINE_TEMPLATE_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned Tq, typename TIn_Type = float>
class CCorrelationEngineTemplate {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;

public:
    virtual const TIn_Type * get_pn() const = 0;

    virtual TIn_Type process(TIn_Type re_in, TIn_Type im_in) = 0;

    virtual ~CCorrelationEngineTemplate() = default;
};

template <unsigned Tq>
class CCorrelationEngineTemplate<Tq, int16_t> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;

public:
#ifndef USE_PN_XOR_TRICK
    virtual const int8_t * get_pn() const = 0;
#else
    virtual const bool * get_pn() const = 0;
#endif
    virtual uint32_t process(int16_t re_in, int16_t im_in) = 0;

    virtual ~CCorrelationEngineTemplate() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORRELATION_ENGINE_TEMPLATE_HPP_
