""" qspannedsequentialadder_fp
"""

# pyright: basic

from aptypes import APComplex, APFixed

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt

from utilities import extract_data, saturate
from qspannedsequentialadder import QSpannedSequentialAdder


class QSpannedSequentialAdderFP:
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

    def __init__(self, q: int, bit_width=16, bit_int=4):
        self._in_width = bit_width
        self._in_int = bit_int
        self._in_quote = bit_width - bit_int
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros((q, 2), dtype=int)

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
        counter = self.__step_counter()
        old_value = self.__fifo[counter].copy()
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
        return value_in - old_value  # pyright: ignore[reportReturnType]


def main():
    GFQ = 64
    N = 60
    RUNS = 5 * 5
    NBR = RUNS * N * GFQ

    IN_W = 8
    IN_I = 6

    var_names = ['data_input_m10dB_w1_q64_N60_0_n30',
                 'iter_fcts_m10dB_w1_q64_N60_0_n30']
    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
    data = dict_data[var_names[0]][0][:NBR]
    qadder_true = dict_data[var_names[1]][0][:NBR]

    DATA = data.astype(np.complex64)

    TMP_MAX = APFixed(0, IN_W, IN_I).max_value
    TMP_MIN = APFixed(0, IN_W, IN_I).min_value
    SATURATED_DATA = saturate(DATA, TMP_MAX, TMP_MIN).astype(
        np.complex64)  # pyright: ignore[reportAttributeAccessIssue]
    del TMP_MIN, TMP_MAX

    FIXED_DATA = np.fromiter((APComplex(z, IN_W, IN_I)
                              for z in SATURATED_DATA), dtype=APComplex)

    qadder = QSpannedSequentialAdder(GFQ)
    qadderfp = QSpannedSequentialAdderFP(GFQ, IN_W, IN_I)

    qadder_resfx: NDArray[np.complex64] = np.array(
        [qadder.process(z.value) for z in FIXED_DATA])
    qadderfp_res: NDArray[np.complex64] = np.array(
        [qadderfp.process(z).value for z in FIXED_DATA])

    plt.figure()
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.title('Real (top) and Imaginary (bottom) parts')
    plt.plot(qadder_true.real, marker='o', label='Pure Float')
    plt.plot(qadder_resfx.real, marker='d', label='Float on Fixed')
    plt.plot(qadderfp_res.real, marker='x', label='Pure Fixed')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(qadder_true.imag, marker='o', label='Pure Float')
    plt.plot(qadder_resfx.imag, marker='d', label='Float on Fixed')
    plt.plot(qadderfp_res.imag, marker='x', label='Pure Fixed')
    plt.legend()

    plt.figure()
    plt.plot(np.abs(qadder_true),  marker='o', label='Pure Float')
    plt.plot(np.abs(qadder_resfx), marker='d', label='Float on Fixed')
    plt.plot(np.abs(qadderfp_res), marker='x', label='Pure Fixed')
    plt.legend()

    plt.show(block=False)

    print(__file__ + ': ok')


if __name__ == "__main__":
    main()
    input('Hit any key to end...')
