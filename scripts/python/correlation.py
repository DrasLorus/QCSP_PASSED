import numpy as np

class StepSub:
    @property
    def q(self):
        return self.__q
    
    @property
    def counter(self):
        return self.__counter

    def rotation(self, n:int):
        return np.exp(2j * self.__omega * n / self.__q)

    def __init__(self, q: int, omega: float):
        self.__omega = omega
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def __step_counter(self):
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1)
        return counter

    def process(self, value_in : np.complex64):
        counter = self.__step_counter()
        new_value = value_in * self.rotation(counter)
        old_value = self.__fifo[counter]
        self.__fifo[counter] = new_value
        return new_value - old_value

class TimeSlidingCorrelator(StepSub):
    @property
    def pn(self):
        return self.__pn

    def __rotatePn(self):
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: np.ndarray, omega: float = 0):
        super().__init__(q, omega)
        self.__pn = pn
        self.__registers = np.zeros(q, dtype=np.complex64)

    def __permute(self, correlations):
        cnt = self.counter
        return np.flip(np.roll(correlations, int(cnt) - 1))


    def __pnCorrelate(self, value_in):
        new_values = value_in * self.pn
        old_correlations = self.__registers
        new_correlations = old_correlations + new_values
        self.__registers = new_correlations
        self.__rotatePn()
        return new_correlations

    def process(self, value_in : np.complex64):
        return self.__pnCorrelate(super().process(value_in))

    def processPermuted(self, value_in : np.complex64):
        return self.__permute(self.process(value_in))

class FftCorrelator:
    @property
    def q(self):
        return self.__q

    @property
    def pn(self):
        return np.fft.ifft(self.__cfpn.conj())

    def rotation(self, n:int):
        return np.exp(2j * self.__omega * n / self.__q).astype(np.complex64)
    
    def __step_counter(self):
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1)
        return counter
    
    def __init__(self, q: int, pn: np.ndarray, omega: float = 0):
        self.__counter = np.uint16(0)
        self.__omega = omega
        self.__q = np.uint16(q)
        self.__cfpn = np.fft.fft(pn).astype(np.complex64).conj()
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def process(self, value_in : np.complex64):
        new_value = value_in * self.rotation(self.__step_counter())
        self.__fifo = np.concatenate((self.__fifo[1:], [new_value]))
        return np.fft.ifft(np.fft.fft(self.__fifo) * self.__cfpn)

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ = 64
    NBR = 10**3 * GFQ
    PN  = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR = 10

    sigma   = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    A = np.array(rng.normal(0., sigma_c, size=NBR) + 1j*rng.normal(0, sigma_c, size=NBR) + np.tile(PN, NBR // GFQ),
                dtype=np.complex64)

    fftCorr = FftCorrelator(GFQ, PN)
    tsCorr  = TimeSlidingCorrelator(GFQ, PN)

    fA = np.array([fftCorr.process(z) for z in A])
    tA = np.array([tsCorr.processPermuted(z) for z in A])

    plt.figure()
    plt.plot(fA[10*GFQ:16*GFQ, 0])
    plt.plot(tA[10*GFQ:16*GFQ, 0])


    plt.figure()
    plt.plot(fA[12*GFQ, :])
    plt.plot(tA[12*GFQ, :])

    plt.figure()
    plt.plot([np.max(np.abs(fA[i-GFQ:i])) for i in range(GFQ,NBR)])
    plt.plot([np.max(np.abs(tA[i-GFQ:i])) for i in range(GFQ,NBR)])

    plt.show()

    print(__file__ + ': ok')