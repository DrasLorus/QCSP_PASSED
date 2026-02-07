""" qspannedsequentialadder
"""

# pyright: basic

import numpy as np
from matplotlib import pyplot as plt

from utilities import extract_data


class QSpannedSequentialAdder:
    """Iterative filter involved in Time sliding-based correlations
    """
    @property
    def q(self) -> np.uint16:
        """value of q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def counter(self) -> np.uint16:
        """getter for the counter

        Returns:
            numpy.uint16: value of the counter
        """
        return self.__counter

    def __init__(self, q: int):
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def reset(self) -> None:
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros(self.q, dtype=np.complex64)

    def __step_counter(self) -> np.uint16:
        """increment and return the value of the counter

        Returns:
            numpy.uint16: value of the stepped counter
        """
        counter = self.__counter
        self.__counter = np.bitwise_and(
            counter + 1, self.__q - 1).astype(np.uint16)
        return counter

    def process(self, value_in: np.complex64) -> np.complex64:
        """process value_in through the iterative filter

        Args:
            value_in (numpy.complex64): new value

        Returns:
            numpy.complex64: resulting value
        """
        counter = self.__step_counter()
        old_value = self.__fifo[counter]
        self.__fifo[counter] = value_in
        return value_in - old_value


def main():
    GFQ = 64
    N = 60
    RUNS = 5 * 5
    NBR = RUNS * N * GFQ

    var_names = ['data_input_m10dB_w1_q64_N60_0_n30',
                 'iter_fcts_m10dB_w1_q64_N60_0_n30']
    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
    data = dict_data[var_names[0]][0][:NBR]
    qadder_true = dict_data[var_names[1]][0][:NBR]

    DATA = data.astype(np.complex64)

    qadder = QSpannedSequentialAdder(GFQ)

    qadder_result = np.array([qadder.process(z) for z in DATA])

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(qadder_result.real, marker='d', label='Python')
    plt.plot(qadder_true.real, marker='x', label='True')
    plt.subplot(2, 1, 2)
    plt.plot(qadder_result.imag, marker='d', label='Python')
    plt.plot(qadder_true.imag, marker='x', label='True')
    plt.legend()

    plt.figure()
    plt.plot(np.abs(qadder_result), marker='d', label='Python')
    plt.plot(np.abs(qadder_true),  marker='x', label='True')
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()

    print(__file__ + ': ok')
