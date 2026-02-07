"""_summary_
"""

# pyright: basic

import numpy as np
from matplotlib import pyplot as plt

from aptypes import APUfixed

from utilities import extract_data
from score_accumulator import ScoreAccumulator


class ScoreAccumulatorFP:
    """score accumulator, fixed-point
    """

    @property
    def q(self) -> np.uint16:
        """getter for self.__q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def p(self) -> np.uint8:
        """getter for self.__p

        Returns:
            int: value of p
        """
        return self.__p

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

    @property
    def bit_width(self) -> int:
        """internal bit_width

        Returns:
            int: the value
        """
        return self.__bit_width

    @property
    def bit_int(self) -> int:
        """internal bit_int

        Returns:
            int: the value
        """
        return self.__bit_int

    def __init__(self, q: int, nsymb: int, bit_width: int, bit_int: int):
        if nsymb > ((2**16 - 1) / q):
            raise ValueError(f"nsymb cannot exceed {2**16 / q - 1}.")
        self.__q = np.uint16(q)
        self.__p = np.uint8(np.log2(q))
        self.__nsymb = np.uint16(nsymb)
        self.__mask = np.uint16(q - 1)
        self.__limit = np.uint16(nsymb * q - 1)
        self.__bit_width = bit_width
        self.__bit_int = bit_int
        self.__fifos = np.array([APUfixed(0, bit_width, bit_int)  # pyright: ignore[reportCallIssue]
                                 for _ in range(q * nsymb)])
        self.__registers = np.array([APUfixed(0, bit_width + int(self.p) + 1, bit_int + int(self.p) + 1)  # pyright: ignore[reportCallIssue]
                                     for _ in range(q)])
        # APUFixed not used to limit performance impact
        self.__counter = np.uint16(0)

    def __step_counter(self) -> tuple[np.uint16, np.uint16]:
        curr_counter = self.__counter
        self.__counter = (curr_counter + 1) * \
            np.uint16(curr_counter != self.__limit)
        return (curr_counter, np.bitwise_and(curr_counter, self.__mask))

    def process(self, value_in: APUfixed) -> APUfixed:
        """process a new max and update internal state

        Args:
            value_in (APUfixed): a new absolute maximum of correlation

        Returns:
            APUfixed: the resulting score
        """
        fifo_counter, register_counter = self.__step_counter()
        old_max = self.__fifos[fifo_counter]
        old_score = self.__registers[register_counter]

        new_score = (old_score + value_in -
                     old_max).truncate(1, False).saturate(1)

        self.__fifos[fifo_counter] = value_in
        self.__registers[register_counter] = new_score

        return new_score


def play(var_names: list[str], data_w: int = 16, data_i: int = 18):
    GFQ = 64
    NFRAME = 60
    RUNS = 10
    NSYMBS = RUNS * (NFRAME * 5)
    NCHIPS = NSYMBS * GFQ

    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names[:2])
    data = dict_data[var_names[0]][0][:NCHIPS]
    true_score = dict_data[var_names[1]][0][:NCHIPS]

    DATA_W = data_w
    DATA_I = data_i

    DATA = data.astype(np.float32)

    TMP_MAX = APUfixed(0, DATA_W, DATA_I).max_value
    TMP_MIN = APUfixed(0, DATA_W, DATA_I).min_value
    SATURATED_DATA = np.clip(DATA, TMP_MIN, TMP_MAX).astype(
        np.float32)  # pyright: ignore[reportAttributeAccessIssue]
    del TMP_MIN, TMP_MAX

    FIXED_DATA = np.fromiter((APUfixed(x, DATA_W, DATA_I)
                              for x in SATURATED_DATA), dtype=APUfixed)

    accufl = ScoreAccumulator(GFQ, NFRAME)
    accufl_scorefp = np.array([accufl.process(value.value)
                              for value in FIXED_DATA])

    accufp = ScoreAccumulatorFP(GFQ, NFRAME, bit_width=DATA_W, bit_int=DATA_I)
    accufp_scorefp = np.array([accufp.process(value) for value in FIXED_DATA])

    fig = plt.figure()
    plt.title(f'Scenario {var_names[2]}')
    plt.plot(true_score, '-o', label='Pure Float')
    plt.plot(accufl_scorefp, '-d', label='Fixed on Float')
    plt.plot(accufp_scorefp, '-x', label='Pure Fixed')
    plt.legend()
    return fig


def main():
    var_names: dict[str, list[str]] = dict()
    var_names['raw_m10dB'] = ['cabs_max_raw_m10dB_w1_q64_N60_0_n30',
                              'score_raw_m10dB_w1_q64_N60_0_n30']
    var_names['l2_m10dB'] = ['cabs_max_l2_m10dB_w1_q64_N60_0_n30',
                             'score_l2_m10dB_w1_q64_N60_0_n30']
    var_names['sqr_raw_m10dB'] = ['cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30',
                                  'score_sqr_raw_m10dB_w1_q64_N60_0_n30']
    var_names['sqr_l2_m10dB'] = ['cabs_max_sqr_l2_m10dB_w1_q64_N60_0_n30',
                                 'score_sqr_l2_m10dB_w1_q64_N60_0_n30']
    var_names['raw_infdB'] = ['cabs_max_raw_infdB_w1_q64_N60_0_n10',
                              'score_raw_infdB_w1_q64_N60_0_n10']
    var_names['l2_infdB'] = ['cabs_max_l2_infdB_w1_q64_N60_0_n10',
                             'score_l2_infdB_w1_q64_N60_0_n10']
    var_names['sqr_raw_infdB'] = ['cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10',
                                  'score_sqr_raw_infdB_w1_q64_N60_0_n10']
    var_names['sqr_l2_infdB'] = ['cabs_max_sqr_l2_infdB_w1_q64_N60_0_n10',
                                 'score_sqr_l2_infdB_w1_q64_N60_0_n10']

    _ = [play(variables + [key]) for key, variables in var_names.items()]
    plt.show(block=False)
    plt.close('all')


if __name__ == '__main__':
    main()
    input('Hit any key to end...')
