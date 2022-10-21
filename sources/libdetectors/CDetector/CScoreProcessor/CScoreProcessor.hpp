#ifndef _QCSP_PASSED_SCORE_PROCESSOR_HPP_
#define _QCSP_PASSED_SCORE_PROCESSOR_HPP_ 1

#include <cstdint>
#include <cstring>
#include <memory>

#include "Miscellanous/misc.hpp"

#include "./CCorrelationEngine/CCorrelationEngine.hpp"
#include "./CNorm/CNorm.hpp"
#include "./CScoreAccumulator/CScoreAccumulator.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type = float, bool normed = true, CorrelationEngineType variant = TIME_SLIDING>
class CScoreProcessor;

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type, CorrelationEngineType variant>
class CScoreProcessor<TFrameSize, Tq, TIn_Type, true, variant> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned N = TFrameSize;

private:
    CCorrelationEngine<Tq, TIn_Type, variant> corr_engine;

    CNorm<q, TIn_Type>                norm_proc;
    CScoreAccumulator<N, q, TIn_Type> score_accumulator;

public:
    const TIn_Type * get_pn() const { return corr_engine.get_pn(); }

    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {
        const TIn_Type norm_unsafe = norm_proc.process(re_in, im_in);
        const TIn_Type norm_value  = norm_unsafe < TIn_Type(1e-8)
                                       ? TIn_Type(1)
                                       : norm_unsafe;

        const TIn_Type cabs_max = std::sqrt(corr_engine.process(re_in, im_in) / norm_value);

        return score_accumulator.process(cabs_max);
    }

    TIn_Type process_sqr(TIn_Type re_in, TIn_Type im_in) {
        const TIn_Type norm_unsafe = norm_proc.process(re_in, im_in);
        const TIn_Type norm_value  = norm_unsafe < TIn_Type(1e-8)
                                       ? TIn_Type(1)
                                       : norm_unsafe;

        const TIn_Type cabs_max = corr_engine.process(re_in, im_in) / norm_value;

        return score_accumulator.process(cabs_max);
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : corr_engine(pn),
          norm_proc(),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type, CorrelationEngineType variant>
class CScoreProcessor<TFrameSize, Tq, TIn_Type, false, variant> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned N = TFrameSize;

private:
    CCorrelationEngine<q, TIn_Type, variant> corr_engine;

    CScoreAccumulator<N, q, TIn_Type> score_accumulator;

public:
    const TIn_Type * get_pn() const { return corr_engine.get_pn(); }

    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {

        const TIn_Type cabs_max = std::sqrt(corr_engine.process(re_in, im_in));

        return score_accumulator.process(cabs_max);
    }

    TIn_Type process_sqr(TIn_Type re_in, TIn_Type im_in) {

        const TIn_Type cabs_max = corr_engine.process(re_in, im_in);

        return score_accumulator.process(cabs_max);
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : corr_engine(pn),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

template <unsigned TFrameSize, unsigned Tq, CorrelationEngineType variant>
class CScoreProcessor<TFrameSize, Tq, int16_t, true, variant> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned p = pow2_log2<q>();
    static constexpr unsigned N = TFrameSize;

private:
    CCorrelationEngine<q, int16_t, variant> corr_engine;

    CNorm<q, int16_t>                norm_proc;
    CScoreAccumulator<N, q, int16_t> score_accumulator;

public:
#ifndef USE_PN_XOR_TRICK
    const int8_t * get_pn() const { return corr_engine.get_pn(); }
#else
    const bool * get_pn() const { return corr_engine.get_pn(); }
#endif

    uint32_t process_sqr(int16_t re_in, int16_t im_in) {

        const uint16_t norm_value = std::max(norm_proc.process_sqr(re_in, im_in), uint16_t(1)); // This max is not useful in a real case scenario

        const uint64_t cabs_max = corr_engine.process(re_in, im_in);

        const uint64_t normed_max = cabs_max / norm_value;

        const uint32_t trunc_mno = normed_max >> (p - 1);

        const uint32_t score = score_accumulator.process(trunc_mno);

        //* To deeply debug, include <cstdio> and uncomment the following
        // fprintf(stderr, "%-8u %-8u %-8u %-8lu %-8lu %-8u %-8u\n", re_in, im_in, norm_value, cabs_max, normed_max, trunc_mno, score);

        return score;
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : corr_engine(pn),
          norm_proc(),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_SCORE_PROCESSOR_HPP_
