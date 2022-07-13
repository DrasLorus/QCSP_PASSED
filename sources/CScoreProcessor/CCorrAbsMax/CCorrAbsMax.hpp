#ifndef _QCSP_PASSED_CORR_ABS_MAX_HPP_
#define _QCSP_PASSED_CORR_ABS_MAX_HPP_ 1

#include <cstdint>
#include <cstring>

#include "Miscellanous/misc.hpp"

namespace QCSP {
namespace StandaloneDetector {

template <unsigned Tq, typename TIn_Type = float>
class CCorrAbsMax {
    static_assert(is_pow2(Tq), "q must be a power of 2.");

public:
    static constexpr unsigned q    = Tq;
    static constexpr unsigned mask = q - 1;

private:
    TIn_Type pn[q];
    TIn_Type re_corr_registers[q];
    TIn_Type im_corr_registers[q];
    TIn_Type abs_corr_registers[q];

    uint32_t counter;

public:
    TIn_Type process(TIn_Type re_in, TIn_Type im_in) {
        const uint32_t curr_counter = counter++;
        for (unsigned u = 0; u < q; u++) {
            const TIn_Type re_pn_u     = pn[(u + curr_counter) & mask];
            const TIn_Type re_corr_i_u = re_corr_registers[u];
            const TIn_Type im_corr_i_u = im_corr_registers[u];

            // local_corr_i[u] = local_buff_i[u] + cpx_pn_u * iterative_factor;
            const TIn_Type re_tmp_corr = re_corr_i_u + re_pn_u * re_in;
            const TIn_Type im_tmp_corr = im_corr_i_u + re_pn_u * im_in;

            re_corr_registers[u] = re_tmp_corr;
            im_corr_registers[u] = im_tmp_corr;

            abs_corr_registers[u] = re_tmp_corr * re_tmp_corr + im_tmp_corr * im_tmp_corr;
        }

        return max_pow2<q>::max(abs_corr_registers);
    }

    template <typename Tpn>
    CCorrAbsMax(Tpn * _pn)
        : counter(0) {
        std::copy(_pn, _pn + q, pn);

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
    static constexpr unsigned mask = q - 1;

private:
    int8_t   pn[q];
    int32_t  re_corr_registers[q];
    int32_t  im_corr_registers[q];
    uint64_t abs_corr_registers[q];

    uint32_t counter;

public:
    uint32_t process(int32_t re_in, int32_t im_in) { // Inputs can be quantified on 17 bits
#ifdef USE_PN_XOR_TRICK
        static constexpr uint32_t transit[2] = {0x0, 0xFFFFFFFF};
#endif
        const uint32_t curr_counter = counter++;
        for (unsigned u = 0; u < q; u++) {

#ifdef USE_PN_XOR_TRICK
            const bool pn_u = pn[(u + curr_counter) & mask];
#else
            const int8_t  pn_u        = pn[(u + curr_counter) & mask];
#endif
            const int32_t re_corr_i_u = re_corr_registers[u]; // Correlations are on 17 + log2q bits
            const int32_t im_corr_i_u = im_corr_registers[u]; // Correlations are on 17 + log2q bits

#ifdef USE_PN_XOR_TRICK
            const uint32_t transit_u = transit[pn_u];

            // local_corr_i[u] = local_buff_i[u] + cpx_pn_u * iterative_factor;
            const int32_t re_tmp_corr = re_corr_i_u + int32_t(re_in ^ transit_u) + pn_u;
            const int32_t im_tmp_corr = im_corr_i_u + int32_t(im_in ^ transit_u) + pn_u;
#else
            const int32_t re_tmp_corr = re_corr_i_u + re_in * pn_u;
            const int32_t im_tmp_corr = im_corr_i_u + im_in * pn_u;
#endif

            re_corr_registers[u] = re_tmp_corr;
            im_corr_registers[u] = im_tmp_corr;
            
            // abs_coor is on 2 * (17 + log2q) + 1 bits, meaning q could be up to 16384 before an overflow
            abs_corr_registers[u] = uint64_t(re_tmp_corr * re_tmp_corr)
                                  + uint64_t(im_tmp_corr * im_tmp_corr); 
        }

        return max_pow2<q>::max(abs_corr_registers);
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
        memset(abs_corr_registers, 0, sizeof(uint64_t) * q);
    }

    virtual ~CCorrAbsMax() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORR_ABS_MAX_HPP_
