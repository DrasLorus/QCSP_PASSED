#ifndef _QCSP_PASSED_SCORE_PROCESSOR_HPP_
#define _QCSP_PASSED_SCORE_PROCESSOR_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

#include "./CCorrAbsMax/CCorrAbsMax.hpp"
#include "./CIterativeAdder/CIterativeAdder.hpp"
#include "./CNorm/CNorm.hpp"
#include "./CScoreAccumulator/CScoreAccumulator.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type = float, bool normed = true>
class CScoreProcessor;

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type>
class CScoreProcessor<TFrameSize, Tq, TIn_Type, true> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned N = TFrameSize;

private:
    CIterativeAdder<q, TIn_Type>      iterative_adder;
    CNorm<q, TIn_Type>                norm_proc;
    CCorrAbsMax<q, TIn_Type>          corr_abs_max;
    CScoreAccumulator<N, q, TIn_Type> score_accumulator;

public:
    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {
        TIn_Type re_it, im_it;

        const TIn_Type norm_unsafe = norm_proc.process(re_in, im_in);
        const TIn_Type norm_value  = norm_unsafe < TIn_Type(1e-8)
                                       ? TIn_Type(1)
                                       : norm_unsafe;

        iterative_adder.process(re_in, im_in, &re_it, &im_it);

        const TIn_Type cabs_max = std::sqrt(corr_abs_max.process(re_it, im_it) / norm_value);

        return score_accumulator.process(cabs_max);
    }

    TIn_Type process_sqr(TIn_Type re_in, TIn_Type im_in) {
        TIn_Type re_it, im_it;

        const TIn_Type norm_unsafe = norm_proc.process(re_in, im_in);
        const TIn_Type norm_value  = norm_unsafe < TIn_Type(1e-8)
                                       ? TIn_Type(1)
                                       : norm_unsafe;

        iterative_adder.process(re_in, im_in, &re_it, &im_it);

        const TIn_Type cabs_max = corr_abs_max.process(re_it, im_it) / norm_value;

        return score_accumulator.process(cabs_max);
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : iterative_adder(),
          norm_proc(),
          corr_abs_max(pn),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

template <unsigned TFrameSize, unsigned Tq, typename TIn_Type>
class CScoreProcessor<TFrameSize, Tq, TIn_Type, false> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned N = TFrameSize;

private:
    CIterativeAdder<q, TIn_Type>      iterative_adder;
    CCorrAbsMax<q, TIn_Type>          corr_abs_max;
    CScoreAccumulator<N, q, TIn_Type> score_accumulator;

public:
    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {
        TIn_Type re_it, im_it;

        iterative_adder.process(re_in, im_in, &re_it, &im_it);

        const TIn_Type cabs_max = std::sqrt(corr_abs_max.process(re_it, im_it));

        return score_accumulator.process(cabs_max);
    }

    TIn_Type process_sqr(TIn_Type re_in, TIn_Type im_in) {
        TIn_Type re_it, im_it;

        iterative_adder.process(re_in, im_in, &re_it, &im_it);

        const TIn_Type cabs_max = corr_abs_max.process(re_it, im_it);

        return score_accumulator.process(cabs_max);
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : iterative_adder(),
          corr_abs_max(pn),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

template <unsigned TFrameSize, unsigned Tq>
class CScoreProcessor<TFrameSize, Tq, int16_t, true> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q = Tq;
    static constexpr unsigned N = TFrameSize;

private:
    CIterativeAdder<q, int16_t>      iterative_adder;
    CNorm<q, int16_t>                norm_proc;
    CCorrAbsMax<q, int16_t>          corr_abs_max;
    CScoreAccumulator<N, q, int16_t> score_accumulator;

public:
    uint64_t process_sqr(int16_t re_in, int16_t im_in) {
        int32_t re_it, im_it;

        const uint32_t norm_value = std::max(norm_proc.process_sqr(re_in, im_in), 1U << 8);

        iterative_adder.process(re_in, im_in, &re_it, &im_it);

        const uint32_t cabs_max = corr_abs_max.process(re_it, im_it);

        const uint32_t normed_max = cabs_max / (norm_value >> 8); // Remove some norm to keep power in the correlation

        return score_accumulator.process(normed_max);
    }

    template <typename Tpn>
    CScoreProcessor(Tpn * pn)
        : iterative_adder(),
          norm_proc(),
          corr_abs_max(pn),
          score_accumulator() {
    }

    virtual ~CScoreProcessor() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_SCORE_PROCESSOR_HPP_
