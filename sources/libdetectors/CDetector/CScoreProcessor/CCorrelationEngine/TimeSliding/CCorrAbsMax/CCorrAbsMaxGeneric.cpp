#include <algorithm>
#include <cstdint>
#include <vector>

#include "./CCorrAbsMaxGeneric.hpp"

float QCSP::StandaloneDetector::CCorrAbsMaxGeneric::process(float re_in, float im_in) {
    for (unsigned u = 0; u < q; u++) {
        const float re_pn_u     = rotating_pn[u];
        const float re_corr_i_u = re_corr_registers[u];
        const float im_corr_i_u = im_corr_registers[u];

        // local_corr_i[u] = local_buff_i[u] + cpx_pn_u * iterative_factor;
        const float re_tmp_corr = re_corr_i_u + re_pn_u * re_in;
        const float im_tmp_corr = im_corr_i_u + re_pn_u * im_in;

        re_corr_registers[u] = re_tmp_corr;
        im_corr_registers[u] = im_tmp_corr;

        abs_corr_registers[u] = re_tmp_corr * re_tmp_corr + im_tmp_corr * im_tmp_corr;
    }

    std::rotate(rotating_pn.begin(), rotating_pn.begin() + 1, rotating_pn.end());

    return *std::max_element(abs_corr_registers.begin(), abs_corr_registers.end());
}

QCSP::StandaloneDetector::CCorrAbsMaxGeneric::CCorrAbsMaxGeneric(const std::vector<float> & _pn)
    : q(_pn.size()),
      mask(q - 1),
      pn(_pn),
      rotating_pn(pn),
      re_corr_registers(q, 0),
      im_corr_registers(q, 0),
      abs_corr_registers(q, 0) {
}
