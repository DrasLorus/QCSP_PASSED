"""exploration of fixed-point division methods
"""

# pyright: basic

from aptypes import APUfixed, APFixed, APComplex

from utilities import saturate
from corrabsmax_fp import TSCorrAbsMaxFP, TSCorrAbsMax
from norm_fp import NormFP, Norm

import numpy as np

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ = int(256)
    RUNS = 50 // 50
    NFRM = 1
    NBR = RUNS * (NFRM + 2) * GFQ
    PN = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR = -10

    sigma = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    data = np.array((rng.normal(0., sigma_c, size=NBR) + 1j*rng.normal(0, sigma_c, size=NBR)
                     + np.tile(np.concatenate((np.zeros(GFQ * NFRM),
                                               np.tile(PN, NFRM),
                                               np.zeros(GFQ * NFRM))), RUNS)),
                    dtype=np.complex64) / 2
    
    IN_W = 8
    IN_I = 2

    data = data * 2**(IN_I - 2) / max(np.max(np.abs(data.real)),
                                                 np.max(np.abs(data.imag)))


    max_temp = APFixed(0, IN_W, IN_I).max_value
    min_temp = APFixed(0, IN_W, IN_I).min_value
    saturated_data = np.array([saturate(a, max_temp, min_temp) for a in data])
    FIXED_DATA = np.array([APComplex(z, IN_W, IN_I) for z in saturated_data], dtype=APComplex)

    del min_temp, max_temp

    ts_corr_fp = TSCorrAbsMaxFP(GFQ, PN, IN_W, IN_I)
    xcam_fp_sat = np.array(
        [ts_corr_fp.process_sqr(z) for z in FIXED_DATA])
    out_width = xcam_fp_sat[0].bit_width
    out_int = xcam_fp_sat[0].bit_int

    ts_corr = TSCorrAbsMax(GFQ, PN)
    xcam_flt_sat = np.array([ts_corr.process_sqr(np.complex64(z.value)) for z in FIXED_DATA])

    norm_calc_fp = NormFP(GFQ, IN_W, IN_I)
    norm_fp = np.array([norm_calc_fp.process(z)
                        for z in FIXED_DATA])
    TMP_ONE = APUfixed(1, norm_fp[0].bit_width, norm_fp[0].bit_int)
    TMP_NUL = APUfixed(0, norm_fp[0].bit_width, norm_fp[0].bit_int)
    norm_fp[norm_fp == TMP_NUL] = TMP_ONE

    norm_calc_flt = Norm(GFQ)
    norm_flt = np.array([norm_calc_flt.process(np.complex64(z.value)) for z in FIXED_DATA])

    inv_norm_flt_fp = np.fromiter(
        map(lambda x: APUfixed(x, out_width, out_int - norm_fp[GFQ].bit_int // 2),
            1. / norm_fp[GFQ:].astype(np.float32)),
        dtype=np.float32)
    inv_norm_fp_fp = np.fromiter(
        map(lambda x: APUfixed(x, out_width + 1, out_int - norm_fp[GFQ].bit_int + 1),
            APUfixed(1., out_width + norm_fp[0].bit_width, out_int).raw // norm_fp[GFQ:].astype(int)),
        dtype=np.float32)

    xcam_flt_div_norm_flt = xcam_flt_sat[GFQ:] * (1. / norm_flt[GFQ:])
    xcam_fp_div_norm_flt = xcam_fp_sat[GFQ:] * (1. / norm_flt[GFQ:])
    xcam_fp_inv_norm_flt_fp = xcam_fp_sat[GFQ:] * inv_norm_flt_fp
    xcam_fp_inv_norm_fp_fp = xcam_fp_sat[GFQ:] * inv_norm_fp_fp

    # Equivalent to map .pad(norm_fp[0].bit_width) to xcam_fp_sat[GFQ:]
    raw_xcamfp_div_norm_fp = np.fromiter(
        map(lambda x: APUfixed(x, out_width, out_int + norm_fp[0].bit_quote),
            (xcam_fp_sat[GFQ:].astype(int)) // norm_fp[GFQ:].astype(int)),
        dtype=float)

    plt.figure()
    plt.title('Result')
    plt.plot(xcam_flt_div_norm_flt, 'b-', label="xcam_flt_div_norm_flt")
    plt.plot(xcam_fp_div_norm_flt, 'g-', label="xcam_fp_div_norm_flt")
    plt.plot(xcam_fp_inv_norm_flt_fp, 'r-', label="xcam_fp_inv_norm_flt_fp")
    plt.plot(xcam_fp_inv_norm_fp_fp, 'y-', label="xcam_fp_inv_norm_fp_fp")
    plt.plot(raw_xcamfp_div_norm_fp, 'c-', label="raw_xcamfp_div_norm_fp")
    plt.legend()

    plt.figure()
    plt.subplot(2, 2, 1)
    plt.title('Correlation Float Sat')
    plt.plot(xcam_flt_sat, 'b-', label="xcam_flt_sat")
    plt.subplot(2, 2, 2)
    plt.title('Norm Float Sat')
    plt.plot(norm_flt, 'g-', label="norm_flt")
    plt.legend()
    plt.subplot(2, 2, 3)
    plt.title('Correlation FP Sat')
    plt.plot(xcam_fp_sat, 'b-', label="xcam_fp_sat")
    plt.subplot(2, 2, 4)
    plt.title('Norm FP Sat')
    plt.plot(norm_fp, 'g-', label="norm_fp")
    plt.legend()

    plt.figure()
    plt.plot(1 / norm_flt[GFQ:], 'b-', label="1 / norm_flt")
    plt.plot(1 / norm_fp[GFQ:].astype(np.float32), 'r-', label="1 / norm_fp")
    plt.plot(inv_norm_flt_fp, 'g-', label="inv_norm_flt_fp")
    plt.plot(inv_norm_fp_fp, 'y-', label="inv_norm_fp_fp")
    plt.legend()

    plt.show()

    print(raw_xcamfp_div_norm_fp[0])
    print(__file__ + ': ok')
