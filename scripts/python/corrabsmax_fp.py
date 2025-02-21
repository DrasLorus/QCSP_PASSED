from aptypes import APComplex, APFixed


from utilities import saturate
from correlation_fp import TimeSlidingCorrelator, TimeSlidingCorrelatorFP

import numpy as np


class TSCorrAbsMax(TimeSlidingCorrelator):
    """extract the absolute max of correlation
    """

    def process(self, value_in: APComplex) -> np.floating:
        """process value_in through the complete correlation process

        Args:
            value_in (APComplex): new value to process

        Returns:
            np.floating: resulting absolute maxima of correlation.
        """
        return np.max(np.abs(super().process(value_in))**2)


class TSCorrAbsMaxFP(TimeSlidingCorrelatorFP):
    """extract the absolute max of correlation, fixed-point version
    """

    def process(self, value_in: APComplex):
        """process value_in through the complete correlation process

        Args:
            value_in (APComplex): values to process

        Returns:
            np.floating: resulting absolute maxima of correlation.
        """
        full_max_values = max(
            tuple(map(lambda z: z.magn(), super().process(value_in))))
        return full_max_values.saturate(self.p + 1).truncate(self.bit_width + 1)


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    GFQ = 64
    RUNS = 10
    NFRAME = 10
    NSYMBS = RUNS * (NFRAME + (NFRAME // 2) * 2)
    NCHIPS = NSYMBS * GFQ
    PN = np.sign(rng.normal(size=(GFQ))).astype(np.float32)
    SNR = -10.

    sigma = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    DUMMY_FRAME = np.concatenate((np.zeros(GFQ * NFRAME // 2),
                                  np.tile(PN, NFRAME),
                                  np.zeros(GFQ * NFRAME // 2)))

    data = np.array(np.sum(rng.normal(0., sigma_c, size=(NCHIPS, 2)) * (1, 1j), axis=1)
                    + np.tile(DUMMY_FRAME, RUNS),
                    dtype=np.complex64) / 2

    IN_W = 7
    IN_I = 4

    tmp_max = APFixed(0, IN_W, IN_I).max_value
    tmp_min = APFixed(0, IN_W, IN_I).min_value
    saturated_data = np.array([saturate(z, tmp_max, tmp_min) for z in data])
    del tmp_min, tmp_max

    ts_corr_fp = TSCorrAbsMaxFP(GFQ, PN, IN_W, IN_I)
    corr_fp = np.array([ts_corr_fp.process(APComplex(z, IN_W, IN_I))
                       for z in saturated_data])

    ts_corr_flt = TSCorrAbsMax(GFQ, PN)
    corr_flt = np.array([ts_corr_flt.process(z) for z in data])

    ts_corr_flt = TSCorrAbsMax(GFQ, PN)
    corr_flt_sat = np.array([ts_corr_flt.process(z) for z in saturated_data])

    plt.figure()
    plt.title('Absolute value of the correlations')
    plt.plot(corr_flt, 'k-o', label="corr_flt")
    plt.plot(corr_flt_sat, 'b-*', label="corr_flt_sat")
    plt.plot(corr_fp, 'g-*', label="corr_fp")
    plt.legend()

    plt.show()

    # Save data for next test
    with open('data.dat', 'wb') as f:
        f.write(corr_fp.astype(np.float32))

    print(__file__ + ': ok')
