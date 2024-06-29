from aptypes import APFixed,APComplex
from correlation import TimeSlidingCorrelator, FftCorrelator
import numpy as np

class StepSubFP:
    @property
    def q(self) -> np.uint16:
        return self.__q

    @property
    def counter(self) -> np.uint16:
        return self.__counter

    @property
    def bitW(self):
        return self._inW

    @property
    def bitI(self):
        return self._inI

    @property
    def bitQ(self):
        return self._inQ

    def __rotation(self, n:int) -> APComplex:
        rotation = np.exp(2j * self.__omega * (2 * n + 1) / (2 * self.q))
        return APComplex(rotation, self._inW, 2)

    def __init__(self, q: int, omega: float, bitW: int = 16, bitI: int = 4):
        self._inW = bitW
        self._inI = bitI
        self._inQ = bitW - bitI
        self.__omega = omega
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros((q, 2), dtype=int)

    def __step_counter(self) -> np.uint16:
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1)
        return counter

    def process(self, value_in: APComplex) -> tuple:
        counter = self.__step_counter()
        # new_value = (value_in * self.__rotation(counter)).truncate(self._inW - 2).truncate(3, False)
        new_value = value_in
        old_value = APComplex(self.__fifo[counter], self.bitW, self.bitI)
        self.__fifo[counter] = (new_value.real.raw, new_value.imag.raw)
        return new_value - old_value

class TimeSlidingCorrelatorFP(StepSubFP):

    @property
    def p(self):
        return self.__p
    
    @property
    def pn(self):
        return self.__pn

    def __rotatePn(self):
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: np.ndarray, omega: float = 0, bitW: int = 16, bitI: int = 4):
        super().__init__(q, omega, bitW, bitI)
        self.__p = int(np.log2(self.q))
        self.__pn = pn.copy()
        self.__registers = tuple(APComplex(0, self.bitW + self.p + 1, self.bitI + self.p + 1)
                                 for _ in range(q))

    def __permuteCorr(self, corr):
        cnt    = self.counter
        npCorr = np.array(corr, dtype=np.complex64)
        return np.flip(np.roll(npCorr, np.int32(cnt - 1)))

    def __pnCorrelate(self, value_in: APComplex):
        new_values = [(value_in * APFixed(p, 2, 2)).truncate(1, False).saturate(1) for p in self.__pn] # FIXME Les problèmes
        old_corlts = self.__registers
        new_corlts = [(x + y).saturate(1) for x,y in zip(new_values, old_corlts)] # FIXME Ajoute un Add stable
        self.__registers = new_corlts
        self.__rotatePn()
        return new_corlts

    def process(self, value_in: APComplex):
        return self.__pnCorrelate(super().process(value_in))
    
    def processPermuted(self, value_in: APComplex):
        return self.__permuteCorr(self.process(value_in))
    
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
    IN_I = 3

    tmpM = APFixed(0, IN_W, IN_I).max_value
    tmpm = APFixed(0, IN_W, IN_I).min_value
    satA = np.array([min(max(a.real, tmpm), tmpM) + 1j*min(max(a.imag, tmpm), tmpM) for a in A])
    del tmpm, tmpM

    fftCorr  = FftCorrelator(GFQ, PN)
    fA = np.array([fftCorr.process(z) for z in A])

    tsCorrFp = TimeSlidingCorrelatorFP(GFQ, PN, 0, IN_W, IN_I)
    pA = np.array([tsCorrFp.processPermuted(APComplex(z, IN_W, IN_I)) for z in satA])

    tsCorr   = TimeSlidingCorrelator(GFQ, PN, 0)
    tA = np.array([tsCorr.processPermuted(z) for z in A])

    tsCorr   = TimeSlidingCorrelator(GFQ, PN, 0)
    tpA = np.array([tsCorr.processPermuted(z) for z in satA])

    plt.figure()
    plt.title(f'Real (top) and Imaginary (bottom) parts of the first point of correlations between {10*GFQ} and {16*GFQ}')
    plt.subplot(2, 1, 1)
    plt.plot(np.real(tA[10*GFQ:16*GFQ, 0]), label="tA")
    plt.plot(np.real(tpA[10*GFQ:16*GFQ, 0]), label="tpA")
    plt.plot(np.real(pA[10*GFQ:16*GFQ, 0]), label="pA")
    plt.subplot(2, 1, 2)
    plt.plot(np.imag(tA[10*GFQ:16*GFQ, 0]), label="tA")
    plt.plot(np.imag(tpA[10*GFQ:16*GFQ, 0]), label="tpA")
    plt.plot(np.imag(pA[10*GFQ:16*GFQ, 0]), label="pA")
    plt.legend()

    plt.figure()
    plt.title(f'Absolute value of the {12*GFQ}th correlation')
    plt.plot(np.abs(tA[12*GFQ, :]), label="tA")
    plt.plot(np.abs(tpA[12*GFQ, :]), label="tpA")
    plt.plot(np.abs(pA[12*GFQ, :]), label="pA")
    plt.legend()

    plt.figure()
    plt.title('Absolute max of correlation')
    plt.plot([np.max(np.abs(tA[i-GFQ:i])) for i in range(GFQ,NBR)], label="tA")
    plt.plot([np.max(np.abs(tpA[i-GFQ:i])) for i in range(GFQ,NBR)], label="tpA")
    plt.plot([np.max(np.abs(pA[i-GFQ:i])) for i in range(GFQ,NBR)], label="pA")
    plt.legend()

    plt.figure()
    plt.title('Absolute max index of correlation')
    plt.plot([np.argmax(np.abs(tA[i-GFQ:i])) for i in range(GFQ,NBR)],  label="tA")
    plt.plot([np.argmax(np.abs(tpA[i-GFQ:i])) for i in range(GFQ,NBR)], label="tpA")
    plt.plot([np.argmax(np.abs(pA[i-GFQ:i])) for i in range(GFQ,NBR)],  label="pA")
    plt.legend()

    plt.show()

    print(__file__ + ': ok')
