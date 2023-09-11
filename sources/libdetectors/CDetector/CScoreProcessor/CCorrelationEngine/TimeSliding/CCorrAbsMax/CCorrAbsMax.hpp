#ifndef _QCSP_PASSED_CORR_ABS_MAX_HPP_
#define _QCSP_PASSED_CORR_ABS_MAX_HPP_ 1

#include <algorithm>
#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

//* DEBUG: If needed ->
// #include <cstdio>

namespace QCSP {
namespace StandaloneDetector {

template <unsigned Tq, typename TIn_Type = float>
class CCorrAbsMax {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q    = Tq;
    static constexpr unsigned mask = q - 1;

private:
    /**
     * @brief reference PN sequence
     *
     */
    TIn_Type pn[q];

    /**
     * @brief used PN sequence
     *
     */
    TIn_Type rotating_pn[q];

    TIn_Type re_corr_registers[q];
    TIn_Type im_corr_registers[q];
    TIn_Type abs_corr_registers[q];

public:
    const TIn_Type * get_pn() const { return pn; }

    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {
        for (unsigned u = 0; u < q; u++) {
            const TIn_Type re_pn_u     = rotating_pn[u];
            const TIn_Type re_corr_i_u = re_corr_registers[u];
            const TIn_Type im_corr_i_u = im_corr_registers[u];

            // local_corr_i[u] = local_buff_i[u] + cpx_pn_u * iterative_factor;
            const TIn_Type re_tmp_corr = re_corr_i_u + re_pn_u * re_in;
            const TIn_Type im_tmp_corr = im_corr_i_u + re_pn_u * im_in;

            re_corr_registers[u] = re_tmp_corr;
            im_corr_registers[u] = im_tmp_corr;

            abs_corr_registers[u] = re_tmp_corr * re_tmp_corr + im_tmp_corr * im_tmp_corr;
        }

        std::rotate(rotating_pn, rotating_pn + 1, rotating_pn + q);

        return max_pow2<q>::max(abs_corr_registers);
    }

    template <typename Tpn>
    CCorrAbsMax(Tpn * _pn) {
        std::copy(_pn, _pn + q, pn);
        std::copy(_pn, _pn + q, rotating_pn);

        memset(re_corr_registers, 0, sizeof(TIn_Type) * q);
        memset(im_corr_registers, 0, sizeof(TIn_Type) * q);
        memset(abs_corr_registers, 0, sizeof(TIn_Type) * q);
    }

    virtual ~CCorrAbsMax() = default;
};

/**
 * @brief CCorrAbsMax Specialization for 16 bit input in the score processor
 *
 * @details This specialization assume a 16-bit quantified input for the full score processor, as delivered by USRPs
 *
 * @tparam q size of the QCSP sequence
 */
template <unsigned Tq>
class CCorrAbsMax<Tq, int16_t> {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q    = Tq;
    static constexpr unsigned p    = pow2_log2<q>();
    static constexpr unsigned mask = q - 1;
    static constexpr unsigned eta  = 16;

private:
    int8_t   pn[q];
    int32_t  re_corr_registers[q];
    int32_t  im_corr_registers[q];
    uint32_t abs_corr_registers[q];

    uint32_t counter;

public:
    const int8_t * get_pn() const { return pn; }

    uint32_t process(int32_t re_in, int32_t im_in) { // Inputs can be quantized on 17 bits
#ifdef USE_PN_XOR_TRICK
        static constexpr uint32_t transit[2] = {0x0, 0xFFFFFFFF};
#endif
        const uint32_t curr_counter = counter++;
        for (unsigned u = 0; u < q; u++) {

#ifdef USE_PN_XOR_TRICK
            const bool pn_u = pn[(u + curr_counter) & mask];
#else
            const int8_t  pn_u = pn[(u + curr_counter) & mask];
#endif
            const int32_t re_corr_i_u = re_corr_registers[u]; // Correlations are on 17 + log2q bits
            const int32_t im_corr_i_u = im_corr_registers[u]; // Correlations are on 17 + log2q bits

#ifdef USE_PN_XOR_TRICK
            const uint32_t transit_u = transit[pn_u];

            // local_corr_i[u] = local_buff_i[u] + cpx_pn_u * iterative_factor;
            const int32_t re_tmp_corr = re_corr_i_u + int32_t(re_in ^ transit_u) + pn_u;
            const int32_t im_tmp_corr = im_corr_i_u + int32_t(im_in ^ transit_u) + pn_u;
#else
            const int32_t re_d = re_in * pn_u; // Can be 17 + 1 sometimes
            const int32_t im_d = im_in * pn_u; // Can be 17 + 1 sometimes

            const int32_t re_sat_d = (re_d <= -int32_t(1 << (17 - 1)) ? -int32_t(1 << (17 - 1)) + 1 : re_d); // It is thus saturated to 17 bits
            const int32_t im_sat_d = (im_d <= -int32_t(1 << (17 - 1)) ? -int32_t(1 << (17 - 1)) + 1 : im_d); // It is thus saturated to 17 bits

            const int32_t re_tmp_corr = re_corr_i_u + re_sat_d;
            const int32_t im_tmp_corr = im_corr_i_u + im_sat_d;
#endif

            re_corr_registers[u] = re_tmp_corr;
            im_corr_registers[u] = im_tmp_corr;

            // abs_coor is on 2 * (17 + log2q) + 1 bits, meaning q could be up to 16384 before an overflow
            const uint64_t full_abs_corr = uint64_t(int64_t(re_tmp_corr) * int64_t(re_tmp_corr))
                                         + uint64_t(int64_t(im_tmp_corr) * int64_t(im_tmp_corr));

            constexpr uint64_t saturation_th = (1LLU << (2 * eta + p + 2)) - 1LLU; // 40 bits for q = 64

            const bool saturate = full_abs_corr > saturation_th;

            const uint64_t sat_abs_corr = full_abs_corr * uint64_t(not saturate)
                                        + saturation_th * uint64_t(saturate);

            const uint32_t trunc_max = uint32_t(sat_abs_corr >> (eta + 1));

            abs_corr_registers[u] = trunc_max;

            // fprintf(stderr, "%-8lu ", abs_corr_registers[u]);
        }

        const uint32_t max_value = max_pow2<q>::max(abs_corr_registers); // 47 bits for q = 64

        //* To deeply debug, include <cstdio> and uncomment the following
        // fprintf(stderr, "%-8u %-8u %-8lu %-8lu %-8u\n", re_in, im_in, max_value);

        return max_value;
    }

    template <typename Tpn>
    CCorrAbsMax(Tpn * _pn)
        : counter(0) {
#ifdef USE_PN_XOR_TRICK
        for (unsigned u = 0; u < Tq; u++) {
            pn[u] = _pn[u] < 0; // -1 -> true, +1 -> false
        }
#else
        std::copy(_pn, _pn + Tq, pn);
#endif

        memset(re_corr_registers, 0, sizeof(int32_t) * q);
        memset(im_corr_registers, 0, sizeof(int32_t) * q);
        memset(abs_corr_registers, 0, sizeof(uint32_t) * q);
    }

    virtual ~CCorrAbsMax() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORR_ABS_MAX_HPP_
