from aptypes import APFixed, APComplex

from utilities import saturate

from correlation import TimeSlidingCorrelator

import numpy as np

class StepSubFP:
    """perform the 'iterative accumulation'

    """
    @property
    def q(self) -> np.uint16:
        """gf cardinal

        Returns:
            np.uint16: value of q
        """
        return self.__q

    @property
    def counter(self) -> np.uint16:
        """internal counter

        Returns:
            np.uint16: the value
        """
        return self.__counter

    @property
    def bit_width(self) -> int:
        """input quantization

        Returns:
            int: the value
        """
        return self._in_width

    @property
    def bit_int(self):
        """input quantization, int-part

        Returns:
            int: the value
        """
        return self._in_int

    @property
    def bit_quote(self):
        """input quantization, quote-part

        Returns:
            int: the value
        """
        return self._in_quote

    def __init__(self, q: int, bit_width = 16, bit_int = 4):
        self._in_width = bit_width
        self._in_int   = bit_int
        self._in_quote = bit_width - bit_int
        self.__q       = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo    = np.zeros((q, 2), dtype=int)

    def __step_counter(self) -> np.uint16:
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1)
        return counter

    def __update_fifo(self, value_in: APComplex) -> APComplex:
        """update the fifo with the new value and get the oldest one

        Args:
            value_in (APComplex): new value

        Returns:
            APComplex: oldest value in the fifo
        """
        counter              = self.__step_counter()
        old_value            = self.__fifo[counter].copy()
        self.__fifo[counter] = (value_in.real.raw, value_in.imag.raw)
        return APComplex(old_value, self.bit_width, self.bit_int)

    def process(self, value_in: APComplex) -> APComplex:
        """process value_in through the acumulative filter

        Args:
            value_in (APComplex): new value in the filter

        Returns:
            APComplex: the newest value minus the oldest one
        """
        old_value = self.__update_fifo(value_in)
        return value_in - old_value

class TimeSlidingCorrelatorFP(StepSubFP):
    """implement a fixed-point time-slidinhg correlator
    """
    @property
    def p(self) -> int:
        """order of GF

        Returns:
            int: value of p
        """
        return self.__p

    @property
    def pn(self) -> np.ndarray[int]:
        """getter for PN

        Returns:
            np.ndarray[int]: a reference to PN
        """
        return self.__pn

    def __rotate_pn(self):
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: np.ndarray[int], bit_width = 16, bit_int = 4):
        super().__init__(q, bit_width, bit_int)
        self.__p = int(np.log2(self.q))
        self.__pn = pn.copy()
        self.__registers = tuple(
            APComplex(0, self.bit_width + self.p + 1, self.bit_int + self.p + 1)
            for _ in range(q))

    def __permute_corr(self, corr: list[APComplex]) -> np.ndarray[np.complex64]:
        """permute the correlation to match fft correlators and cast to np.complex64

        Args:
            corr (APComplex): correlation to permute

        Returns:
            np.ndarray[np.complex64]: permuted np.complex64-converted correlation
        """
        cnt    = self.counter
        np_corr = np.array(corr, dtype=np.complex64)
        return np.flip(np.roll(np_corr, np.int32(cnt - 1)))

    def __pn_correlate(self, value_in: APComplex) -> list[APComplex]:
        """ts-correlate value in with PN

        Args:
            value_in (APComplex): input value

        Returns:
            list[APComplex]: the new correlation vector
        """
        new_values = [(value_in * APFixed(p, 2, 2)).truncate(1, False).saturate(1)
            for p in self.__pn]
        old_corlts = self.__registers
        new_corlts = [(x + y).saturate(1) for x,y in zip(new_values, old_corlts)]
        self.__registers = new_corlts
        self.__rotate_pn()
        return new_corlts

    def process(self, value_in: APComplex) -> list[APComplex]:
        """process value in through the fixed-point TS-based correlator

        Args:
            value_in (APComplex): new value

        Returns:
            list[APComplex]: correlation vector
        """
        return self.__pn_correlate(super().process(value_in))

    def process_permuted(self, value_in: APComplex):
        """retrieve the permuted correlation, matching FFT-based ones

        Args:
            value_in (APComplex): new value

        Returns:
            list[APComplex]: permuted correlation vector
        """
        return self.__permute_corr(self.process(value_in))

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    GFQ = 64
    NBR = 100 * GFQ
    PN  = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR = -10

    sigma   = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    data = np.array(np.sum(rng.normal(0., sigma_c, size=(NBR, 2)) * (1, 1j), axis=1)
                 + np.tile(PN, NBR // GFQ),
                dtype=np.complex64) / 2

    IN_W = 7
    IN_I = 3

    tmp_max        = APFixed(0, IN_W, IN_I).max_value
    tmp_min        = APFixed(0, IN_W, IN_I).min_value
    saturated_data = np.array([saturate(z, tmp_max, tmp_min) for z in data])
    del tmp_min, tmp_max

    # fft_corr_proc  = FftCorrelator(GFQ, PN)
    # fft_corr       = np.array([fft_corr_proc.process(z) for z in data])

    ts_corr_proc_fp = TimeSlidingCorrelatorFP(GFQ, PN, IN_W, IN_I)
    ts_corr_fp      = np.array([ts_corr_proc_fp.process_permuted(APComplex(z, IN_W, IN_I))
                                for z in saturated_data])

    ts_corr_proc_flt = TimeSlidingCorrelator(GFQ, PN)
    ts_corr_flt      = np.array([ts_corr_proc_flt.process_permuted(z) for z in data])

    ts_corr_proc_flt = TimeSlidingCorrelator(GFQ, PN)
    ts_corr_flt_sat  = np.array([ts_corr_proc_flt.process_permuted(z) for z in saturated_data])

    plt.figure()
    plt.title('Real (top) and Imaginary (bottom) parts of the first point of correlations between '
              + f'{10*GFQ} and {16*GFQ}')
    plt.subplot(2, 1, 1)
    plt.plot(np.real(ts_corr_flt[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_flt")
    plt.plot(np.real(ts_corr_flt_sat[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_flt_sat")
    plt.plot(np.real(ts_corr_fp[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_fp")
    plt.subplot(2, 1, 2)
    plt.plot(np.imag(ts_corr_flt[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_flt")
    plt.plot(np.imag(ts_corr_flt_sat[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_flt_sat")
    plt.plot(np.imag(ts_corr_fp[10*GFQ:16*GFQ, 0]),
        label = "ts_corr_fp")
    plt.legend()

    plt.figure()
    plt.title(f'Absolute value of the {12*GFQ}th correlation')
    plt.plot(np.abs(ts_corr_flt[12*GFQ, :]),
        label = "ts_corr_flt", marker='*')
    plt.plot(np.abs(ts_corr_flt_sat[12*GFQ, :]),
        label = "ts_corr_flt_sat", marker='o')
    plt.plot(np.abs(ts_corr_fp[12*GFQ, :]),
        label = "ts_corr_fp", marker='x')
    plt.legend()

    plt.figure()
    plt.title('Absolute max of correlation')
    plt.plot([np.max(np.abs(ts_corr_flt[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_flt")
    plt.plot([np.max(np.abs(ts_corr_flt_sat[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_flt_sat")
    plt.plot([np.max(np.abs(ts_corr_fp[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_fp")
    plt.legend()

    plt.figure()
    plt.title('Absolute max index of correlation')
    plt.plot([np.argmax(np.abs(ts_corr_flt[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_flt")
    plt.plot([np.argmax(np.abs(ts_corr_flt_sat[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_flt_sat")
    plt.plot([np.argmax(np.abs(ts_corr_fp[i-GFQ:i])) for i in range(GFQ,NBR)],
        label = "ts_corr_fp")
    plt.legend()

    plt.show()

    print(__file__ + ': ok')
