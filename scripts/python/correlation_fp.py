"""Comparison of correlation method approach
"""

# pyright: basic

from argparse import ArgumentParser

from matplotlib import pyplot as plt
from numpy.typing import NDArray
import numpy as np

from aptypes import APFixed, APComplex

from utilities import saturate, generate_data, extract_data

from correlation import TimeSlidingCorrelator

from qspannedsequentialadder_fp import QSpannedSequentialAdderFP


class TimeSlidingCorrelatorFP:
    """implement a fixed-point time-sliding correlator
    """
    @property
    def p(self) -> np.uint8:
        """order of GF

        Returns:
            int: value of p
        """
        return self.__p

    @property
    def q(self) -> np.uint16:
        """order of GF

        Returns:
            int: value of p
        """
        return self.__adder.q

    @property
    def bit_width(self) -> int:
        """input quantization

        Returns:
            int: the value
        """
        return self.__adder.bit_width

    @property
    def bit_int(self) -> int:
        """input quantization, int-part

        Returns:
            int: the value
        """
        return self.__adder.bit_int

    @property
    def pn(self) -> NDArray[np.int8]:
        """getter for PN

        Returns:
            NDArray[int]: a reference to PN
        """
        return self.__pn

    def __rotate_pn(self):
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: NDArray[np.integer | np.floating], bit_width: int = 16, bit_int: int = 6):
        self.__adder = QSpannedSequentialAdderFP(q, bit_width, bit_int)
        self.__p: np.uint8 = np.log2(self.q).astype(np.uint8)
        self.__pn = pn.astype(np.int8)
        self.__registers = list(
            APComplex(0, int(self.bit_width + self.p + 1),
                      int(self.bit_int + self.p + 1))
            for _ in range(q))

    def __permute_corr(self, corr: list[APComplex]) -> NDArray[np.complex64]:
        """permute the correlation to match fft correlators and cast to np.complex64

        Args:
            corr (APComplex): correlation to permute

        Returns:
            NDArray[np.complex64]: permuted np.complex64-converted correlation
        """
        cnt = self.__adder.counter
        np_corr = np.array(corr, dtype=np.complex64)
        return np.flip(np.roll(np_corr, np.int32(cnt - 1)))

    def __pn_correlate(self, value_in: APComplex) -> list[APComplex]:
        """ts-correlate value in with PN

        Args:
            value_in (APComplex): input value

        Returns:
            list[APComplex]: the new correlation vector
        """
        new_values: list[APComplex] = [(value_in * APFixed(p, 2, 2)).truncate(1, False).saturate(1)
                                       for p in self.__pn]
        old_corlts: list[APComplex] = self.__registers
        new_corlts: list[APComplex] = [
            (x + y).saturate(1) for x, y in zip(new_values, old_corlts)]
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
        return self.__pn_correlate(self.__adder.process(value_in))

    def process_permuted(self, value_in: APComplex) -> NDArray[np.complex64]:
        """retrieve the permuted correlation, matching FFT-based ones

        Args:
            value_in (APComplex): new value

        Returns:
            list[APComplex]: permuted correlation vector
        """
        return self.__permute_corr(self.process(value_in))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--data', type=str,
                        default='generate', choices=['generate', 'extract'])
    args = parser.parse_args()

    GFQ = 64
    NFM = 60
    RUNS = NFM * 3
    NBR = RUNS * GFQ
    SNR = -10.
    SEED = 0

    DATA_METHOD: str = args.data

    EXTRACTED = DATA_METHOD == 'extract'
    GENERATED = DATA_METHOD == 'generate'

    if GENERATED:
        pn = np.sign(np.random.randn(GFQ))
        data = generate_data(GFQ, RUNS, SNR, pn, SEED)
    elif EXTRACTED:
        pn = extract_data('data/parameters_20210903.mat',
                          ['PN64'])['PN64'].reshape((64))
        var_names = ['data_input_m10dB_w1_q64_N60_0_n30']
        dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
        data = dict_data['data_input_m10dB_w1_q64_N60_0_n30'][0][:NBR]
    DATA = data.astype(np.complex64)
    PN = pn.astype(np.float32)

    IN_W = 8
    IN_I = 6

    tmp_max = APFixed(0, IN_W, IN_I).max_value
    tmp_min = APFixed(0, IN_W, IN_I).min_value
    saturated_data = saturate(data, tmp_max, tmp_min).astype(
        np.complex64)  # pyright: ignore[reportAttributeAccessIssue]
    del tmp_min, tmp_max

    fixed_data = np.fromiter(
        (APComplex(z, IN_W, IN_I).value for z in saturated_data), dtype=np.complex64)

    ts_corr_proc_fp = TimeSlidingCorrelatorFP(GFQ, PN, IN_W, IN_I)
    ts_corr_fp = np.array([ts_corr_proc_fp.process_permuted(APComplex(z, IN_W, IN_I))
                           for z in saturated_data])

    ts_corr_proc_flt = TimeSlidingCorrelator(GFQ, PN)
    ts_corr_flt = np.array(
        [ts_corr_proc_flt.process_permuted(z) for z in data])

    ts_corr_proc_flt.reset()
    ts_corr_flt_sat = np.array(
        [ts_corr_proc_flt.process_permuted(z) for z in fixed_data])

    plt.figure()
    plt.title('Real (top) and Imaginary (bottom) parts of the first point of correlations between '
              + f'{10*GFQ} and {16*GFQ}')
    plt.subplot(2, 1, 1)
    plt.plot(np.real(ts_corr_flt[10*GFQ:16*GFQ, 0]), label="ts_corr_flt")
    plt.plot(
        np.real(ts_corr_flt_sat[10*GFQ:16*GFQ, 0]), label="ts_corr_flt_sat")
    plt.plot(np.real(ts_corr_fp[10*GFQ:16*GFQ, 0]), label="ts_corr_fp")
    plt.subplot(2, 1, 2)
    plt.plot(np.imag(ts_corr_flt[10*GFQ:16*GFQ, 0]), label="ts_corr_flt")
    plt.plot(
        np.imag(ts_corr_flt_sat[10*GFQ:16*GFQ, 0]), label="ts_corr_flt_sat")
    plt.plot(np.imag(ts_corr_fp[10*GFQ:16*GFQ, 0]), label="ts_corr_fp")
    plt.legend()

    plt.figure()
    plt.title(f'Absolute value of the {12*GFQ}th correlation')
    plt.plot(np.abs(ts_corr_flt[12*GFQ, :]), label="ts_corr_flt", marker='*')
    plt.plot(np.abs(ts_corr_flt_sat[12*GFQ, :]),
             label="ts_corr_flt_sat", marker='o')
    plt.plot(np.abs(ts_corr_fp[12*GFQ, :]), label="ts_corr_fp", marker='x')
    plt.legend()

    plt.figure()
    plt.title('Absolute max of correlation')
    plt.plot(np.max(np.abs(ts_corr_flt), axis=1), label="Pure Float")
    plt.plot(np.max(np.abs(ts_corr_flt_sat), axis=1), label="Float on Fixed")
    plt.plot(np.max(np.abs(ts_corr_fp), axis=1), label="Pure Fixed")
    plt.legend()

    plt.figure()
    plt.title('Absolute max index of correlation')
    plt.plot(np.argmax(np.abs(ts_corr_flt), axis=1), label="ts_corr_flt")
    plt.plot(np.argmax(np.abs(ts_corr_flt_sat), axis=1),
             label="ts_corr_flt_sat")
    plt.plot(np.argmax(np.abs(ts_corr_fp), axis=1), label="ts_corr_fp")
    plt.legend()

    plt.show()

    print(__file__ + ': ok')
