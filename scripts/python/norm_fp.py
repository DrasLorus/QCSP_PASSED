from aptypes import APFixed, APUfixed, APComplex
import numpy as np

class Norm:
    @property
    def q(self):
        return self.__q

    def __init__(self, q : int):
        self.__q = q
        self._counter = np.uint8(0)
        self._fifo = np.concatenate((np.zeros(q -1 , dtype=np.float32), np.ones(1, dtype=np.float32)))
        self._register = np.float32(1)

    def process(self, value_in : np.complex64):
        new_value = np.real(value_in.conj() * value_in)
        counter = self._counter
        old_value = self._fifo[counter]
        self._fifo[counter] = new_value
        self._counter = np.bitwise_and(counter + 1, self.q - 1)
        old_norm = self._register
        self._register = old_norm + new_value - old_value
        return self._register

class NormFP:
    def __init__(self, q : int, bit_width : int, bit_int : int):
        self._inW = bit_width
        self._inI = bit_int
        self._inQ = bit_width - bit_int
        self._q = q
        self._p = int(np.log2(q))
        self._counter = np.uint8(q - 1)
        self._fifo_qtf = (2 * bit_width + 1 - self._inQ, 2 * bit_int + 1)
        self._fifo = np.concatenate((np.zeros(q - 1), np.ones(1))) # bit_width * 2 + 1, bit_int * 2 + 1
        self._reg_qtf = (2 * bit_width + 2*(self._p + 1) - self._inQ, 2 * bit_int + 2*(self._p + 1))
        self._register = APUfixed(1., *self._reg_qtf)

    def __step_counter(self):
        counter = np.bitwise_and(self._counter + 1, self._q - 1, dtype=np.uint8)
        self._counter = counter
        return counter

    def __step_register(self, new_value):
        counter = self.__step_counter()
        old_value = APUfixed(self._fifo[counter], *self._fifo_qtf)
        self._fifo[counter] = new_value.value
        old_norm = self._register
        self._register = ((old_norm + new_value) - old_value).truncate(1, lsb=False).saturate(1) # Stable
        return self._register

    def process(self, value_in: APComplex):
        new_value = value_in.magn().truncate(self._inQ)
        return self.__step_register(new_value).truncate(self._p + 1)

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ  = int(64)
    RUNS = 2000
    NFRM = 1
    NBR  = RUNS * (NFRM + 2) * GFQ
    PN   = np.sign(np.random.randn(GFQ)).astype(np.float32)

    snr     = -10
    sigma   = np.sqrt(10**(-snr / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    A = np.array((rng.normal(0., sigma_c, size=NBR) + 1j*rng.normal(0, sigma_c, size=NBR) 
                  + np.tile(np.concatenate((np.zeros(64), np.tile(PN, NFRM), np.zeros(64))), RUNS)),
                dtype=np.complex64) / 2

    # PA = np.concatenate((np.zeros(GFQ - 1), A))
    # N  = np.array([np.sum(np.real(PA[i : i + GFQ].conj() * PA[i : i + GFQ]))
    #             for i,_ in enumerate(PA[:-GFQ+1])])

    norm = Norm(GFQ)
    NIT  = np.array([norm.process(x) for x in A])

    IN_W = 8
    IN_I = 4

    tmpM = APFixed(0, IN_W, IN_I).max_value
    tmpm = APFixed(0, IN_W, IN_I).min_value
    SatA = np.array([min(max(z.real, tmpm), tmpM) + 1j*min(max(z.imag, tmpm), tmpM) for z in A])
    del tmpm, tmpM

    norm   = Norm(GFQ)
    normfp = NormFP(GFQ, IN_W, IN_I)

    NITP  = np.array([norm.process(z) for z in SatA])
    NITFP = np.array([normfp.process(APComplex(z, IN_W, IN_I)) for z in SatA])

    plt.figure()
    plt.plot(NIT, 'k:', label="Float")
    plt.plot(NITP, 'b-x', label="Float Sat.")
    plt.plot(NITFP, 'g-x', label="FP Sat.")
    plt.legend()
    plt.show()

    print(__file__ + ': ok')