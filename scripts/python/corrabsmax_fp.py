from aptypes import APComplex, APFixed
from correlation_fp import *

import numpy as np

class TSCorrAbsMax(TimeSlidingCorrelator):
    def process(self, input: np.ndarray):
        return np.max(np.abs(super().process(input))**2)

class TSCorrAbsMaxFP(TimeSlidingCorrelatorFP):
    def process(self, input: APComplex):
        return max(tuple(map(lambda z: z.magn(), super().process(input)))).saturate(self.p + 1).truncate(self.bitW + 1)

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ = 64
    NBR = 100 * GFQ
    PN  = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR = -10

    sigma   = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    A = np.array(rng.normal(0., sigma_c, size=NBR) + 1j*rng.normal(0, sigma_c, size=NBR) + np.tile(PN, NBR // GFQ),
                dtype=np.complex64) / 2

    IN_W = 7
    IN_I = 4

    tmpM = APFixed(0, IN_W, IN_I).maxValue
    tmpm = APFixed(0, IN_W, IN_I).minValue
    satA = np.array([min(max(a.real, tmpm), tmpM) + 1j*min(max(a.imag, tmpm), tmpM) for a in A])
    del tmpm, tmpM

    tsCorrFp = TSCorrAbsMaxFP(GFQ, PN, 0, IN_W, IN_I)
    pA = np.array([tsCorrFp.process(APComplex(z, IN_W, IN_I)) for z in satA])

    tsCorr   = TSCorrAbsMax(GFQ, PN, 0)
    tA = np.array([tsCorr.process(z) for z in A])

    tsCorr   = TSCorrAbsMax(GFQ, PN, 0)
    tpA = np.array([tsCorr.process(z) for z in satA])

    plt.figure()
    plt.title(f'Absolute value of the {12*GFQ}th correlation')
    plt.plot(tA, 'k-o', label="tA")
    plt.plot(tpA, 'b-*', label="tpA")
    plt.plot(pA, 'g-*', label="pA")
    plt.legend()

    plt.show()

    print(__file__ + ': ok')
