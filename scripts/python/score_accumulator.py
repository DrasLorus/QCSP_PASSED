"""_summary_
"""

# pyright: basic

import numpy as np
from matplotlib import pyplot as plt

from utilities import extract_data


class ScoreAccumulator:
    """score accumulator, floating-point
    """

    @property
    def q(self) -> np.uint16:
        """getter for self.__q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def nsymb(self) -> np.uint16:
        """getter for self.__nsymb

        Returns:
            int: value of nsymb
        """
        return self.__nsymb

    @property
    def counter(self) -> np.uint16:
        """internal counter

        Returns:
            np.uint16: the value
        """
        return self.__counter

    def __init__(self, q: int, nsymb: int):
        if nsymb > ((2**16 - 1) / q):
            raise ValueError(f"nsymb cannot exceed {2**16 / q - 1}.")
        self.__q = np.uint16(q)
        self.__nsymb = np.uint16(nsymb)
        self.__mask = np.uint16(q - 1)
        self.__limit = np.uint16(nsymb * q - 1)
        self.__fifos = np.zeros(self.nsymb * self.q, dtype=np.float32)
        self.__registers = np.zeros(self.q, dtype=np.float32)
        self.__counter = np.uint16(0)

    def __step_counter(self) -> tuple[np.uint16, np.uint16]:
        curr_counter = self.__counter
        self.__counter = (curr_counter + 1) * \
            np.uint16(curr_counter != self.__limit)
        return (curr_counter, np.bitwise_and(curr_counter, self.__mask))

    def process(self, value_in: float) -> float:
        """process a new max and update internal state

        Args:
            value_in (float): a new absolute maximum of correlation

        Returns:
            float: the resulting score
        """
        fifo_counter, register_counter = self.__step_counter()
        old_max = self.__fifos[fifo_counter].copy()
        old_score = self.__registers[register_counter].copy()

        new_score = old_score + value_in - old_max

        self.__fifos[fifo_counter] = value_in
        self.__registers[register_counter] = new_score

        return new_score


def main():
    GFQ = 64
    NFRAME = 60
    RUNS = 2
    NSYMBS = RUNS * (NFRAME * 5)
    NCHIPS = NSYMBS * GFQ

    var_names = ['cabs_max_raw_m10dB_w1_q64_N60_0_n30',
                 'score_raw_m10dB_w1_q64_N60_0_n30']
    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
    data = dict_data[var_names[0]][0][:NCHIPS]
    true_result = dict_data[var_names[1]][0][:NCHIPS]

    plt.figure()
    plt.plot(data, label='Data (FP)')

    accu_flt = ScoreAccumulator(GFQ, NFRAME)
    score_flt = np.array([accu_flt.process(value) for value in data])

    plt.figure()
    plt.plot(score_flt, '-o', label='Float')
    plt.plot(true_result, '-x', label='True')
    plt.show()


if __name__ == '__main__':
    main()
