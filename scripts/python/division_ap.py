from aptypes import APUfixed, APFixed, APComplex
from corrabsmax_fp import TSCorrAbsMaxFP, TSCorrAbsMax
from norm_fp import NormFP, Norm

import numpy as np

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ  = int(256)
    RUNS = 5
    NFRM = 1
    NBR  = RUNS * (NFRM + 2) * GFQ
    PN   = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR  = -10

    sigma   = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    A = np.array((rng.normal(0., sigma_c, size=NBR) + 1j*rng.normal(0, sigma_c, size=NBR) 
                  + np.tile(np.concatenate((np.zeros(GFQ * NFRM), np.tile(PN, NFRM), np.zeros(GFQ * NFRM))), RUNS)),
                dtype=np.complex64) / 2

    IN_W = 8
    IN_I = 4

    tmpM = APFixed(0, IN_W, IN_I).maxValue
    tmpm = APFixed(0, IN_W, IN_I).minValue
    satA = np.array([min(max(a.real, tmpm), tmpM) + 1j*min(max(a.imag, tmpm), tmpM) for a in A])
    del tmpm, tmpM

    tsCorrFp = TSCorrAbsMaxFP(GFQ, PN, 0, IN_W, IN_I)
    pA       = np.array([tsCorrFp.process(APComplex(z, IN_W, IN_I)) for z in satA])
    outW = pA[0].bitW
    outI = pA[0].bitI

    tsCorr = TSCorrAbsMax(GFQ, PN, 0)
    tpA    = np.array([tsCorr.process(z) for z in satA])

    normFp = NormFP(GFQ, IN_W, IN_I)
    pB     = np.array([normFp.process(APComplex(z, IN_W, IN_I)) for z in satA])

    norm = Norm(GFQ)
    tpB  = np.array([norm.process(z) for z in satA])

    pBIF = np.fromiter(map(lambda x: APUfixed(x, outW, outI - pB[GFQ].bitI // 2), 1. / pB[GFQ:].astype(np.float32)), dtype=np.float32)
    ivpB = np.fromiter(map(lambda x: APUfixed(x, outW, outI - pB[GFQ].bitI), APUfixed(1., outW + pB[0].bitW, outI).raw // pB[GFQ:].astype(int)), dtype=np.float32)
    
    tpAtpB = tpA[GFQ:] / tpB[GFQ:]
    pAtpB  = pA[GFQ:] * (1. / tpB[GFQ:])
    pApBIF  = pA[GFQ:] * pBIF
    pAivpB  = pA[GFQ:] * ivpB

    pApB = np.fromiter(map(lambda x: APUfixed(x, outW, outI + pB[0].bitQ),
                           (pA[GFQ:].astype(int)) // pB[GFQ:].astype(int)), # Equivalent to map .pad(pB[0].bitW) to pA[GFQ:]
                           dtype=float)

    plt.figure()
    plt.title(f'Result')
    plt.plot(tpAtpB, 'b-', label="tpAtpB")
    plt.plot(pAtpB, 'g-', label="pAtpB")
    plt.plot(pApBIF, 'r-', label="pApBIF")
    plt.plot(pAivpB, 'y-', label="pAivpB")
    plt.plot(pApB, 'c-', label="pApB")
    plt.legend()

    plt.figure()
    plt.subplot(2,2,1)
    plt.title(f'Correlation Float Sat')
    plt.plot(tpA, 'b-', label="tpA")
    plt.subplot(2,2,2)
    plt.title(f'Norm Float Sat')
    plt.plot(tpB, 'g-', label="tpB")
    plt.legend()
    plt.subplot(2,2,3)
    plt.title(f'Correlation FP Sat')
    plt.plot(pA, 'b-', label="pA")
    plt.subplot(2,2,4)
    plt.title(f'Norm FP Sat')
    plt.plot(pB, 'g-', label="pB")
    plt.legend()

    plt.figure()
    plt.plot(1 / tpB[GFQ:], 'b-', label="1 / tpB")
    plt.plot(1 / pB[GFQ:].astype(np.float32), 'r-', label="1 / pB")
    plt.plot(pBIF, 'g-', label="pBIF")
    plt.plot(ivpB, 'y-', label="ivpB")
    plt.legend()

    plt.show()

    print(__file__ + ': ok')
