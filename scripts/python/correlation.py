"""Comparison of correlation method approach
"""

import numpy as np

class StepSub:
    """Iterative filter involved in Time sliding-based correlations
    """
    @property
    def q(self) -> int:
        """value of q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def counter(self) -> np.uint16:
        """getter for the counter

        Returns:
            np.uint16: value of the counter
        """
        return self.__counter

    def __init__(self, q: int):
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def __step_counter(self) -> np.uint16:
        """increment and return the value of the counter

        Returns:
            np.uint16: value of the stepped counter
        """
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1).astype(np.uint16)
        return counter

    def process(self, value_in : np.complex64) -> np.complex64:
        """process value_in through the iterative filter

        Args:
            value_in (np.complex64): new value

        Returns:
            np.complex64: resulting value
        """
        counter = self.__step_counter()
        old_value = self.__fifo[counter]
        self.__fifo[counter] = value_in
        return value_in - old_value

class TimeSlidingCorrelator(StepSub):
    """Time sliding-based correlator, floating-point implementation 
    """
    @property
    def pn(self) -> np.ndarray[np.float32]:
        """the PN sequence

        Returns:
            np.ndarray[np.float32]: the PN sequence
        """
        return self.__pn

    def __rotate_pn(self) -> None:
        """rotate the PN sequence in place
        """
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: np.ndarray):
        super().__init__(q)
        self.__pn = pn
        self.__registers = np.zeros(q, dtype=np.complex64)

    def __permute(self, correlations: np.ndarray[np.complex64]):
        """perform the permutation

        Args:
            correlations (np.ndarray[np.complex64]): correlation array

        Returns:
            np.ndarray[np.complex64]: permuted correlation array
        """
        cnt = self.counter
        return np.flip(np.roll(correlations, int(cnt) - 1))


    def __pn_correlate(self, value_in: np.complex64):
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        new_values       = value_in * self.pn
        old_correlations = self.__registers
        new_correlations = old_correlations + new_values
        self.__registers = new_correlations
        self.__rotate_pn()
        return new_correlations

    def process(self, value_in: np.complex64):
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        return self.__pn_correlate(super().process(value_in))

    def process_permuted(self, value_in : np.complex64) -> np.ndarray[np.complex64]:
        """retrieve the permuted correlation, matching FFT-based ones

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: permuted correlation vector
        """
        return self.__permute(self.process(value_in))

class FftCorrelator:
    """FFT-based correlator, floating-point implementation
    """
    @property
    def q(self) -> np.uint16:
        """getter for q

        Returns:
            np.uint16: the value of q
        """
        return self.__q

    @property
    def pn(self) -> np.ndarray[np.float32]:
        """getter for the PN sequence

        Returns:
            np.ndarray[np.float32]: the sequence
        """
        return np.fft.ifft(self.__cfpn.conj()).real.astype(np.float32)

    def __init__(self, q: int, pn: np.ndarray[np.float32]):
        self.__q = np.uint16(q)
        self.__cfpn = np.fft.fft(pn).astype(np.complex64).conj()
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def process(self, value_in : np.complex64) -> np.ndarray[np.complex64]:
        """perform the correlation using FFT method

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        self.__fifo = np.concatenate((self.__fifo[1:], [value_in]))
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

    data = np.array(np.sum(rng.normal(0.,sigma_c, size=(NBR, 2)) * (1, 1j), axis=1)
                 + np.tile(PN, NBR // GFQ),
                dtype=np.complex64)

    fft_corr_proc = FftCorrelator(GFQ, PN)
    ts_corr_proc  = TimeSlidingCorrelator(GFQ, PN)

    fft_corr = np.array([fft_corr_proc.process(z) for z in data])
    ts_corr  = np.array([ts_corr_proc.process_permuted(z) for z in data])

    plt.figure()
    plt.plot(fft_corr[10*GFQ:16*GFQ, 0], label = 'FFT')
    plt.plot(ts_corr[10*GFQ:16*GFQ, 0], label = 'TS')


    plt.figure()
    plt.plot(fft_corr[12*GFQ, :], label = 'FFT')
    plt.plot(ts_corr[12*GFQ, :], label = 'TS')

    plt.figure()
    plt.plot([np.max(np.abs(fft_corr[i-GFQ:i])) for i in range(GFQ,NBR)], label = 'FFT')
    plt.plot([np.max(np.abs(ts_corr[i-GFQ:i])) for i in range(GFQ,NBR)], label = 'TS')

    plt.show()

    print(__file__ + ': ok')
