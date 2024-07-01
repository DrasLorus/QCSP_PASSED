from utilities import saturate

from aptypes import APUfixed

import numpy as np

class ScoreAccumulator:
    """score accumulator, floating-point
    """

    @property
    def q(self) -> int:
        """getter for self.__q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def nsymb(self) -> int:
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
        self.__q         = q
        self.__nsymb     = nsymb
        self.__mask      = q - 1
        self.__limit     = nsymb * q - 1
        self.__fifos     = np.zeros(self.nsymb * self.q, dtype=np.float32)
        self.__registers = np.zeros(self.q, dtype=np.float32)
        self.__counter   = np.uint16(0)

    def __step_counter(self) -> np.uint16:
        curr_counter = self.__counter
        self.__counter = (curr_counter + 1) * np.uint16(curr_counter != self.__limit)
        return (curr_counter, np.bitwise_and(curr_counter, self.__mask))

    def process(self, value_in: float) -> float:
        """process a new max and update internal state

        Args:
            value_in (float): a new absolute maximum of correlation

        Returns:
            float: the resulting score
        """
        fifo_counter, register_counter = self.__step_counter()
        old_max   = self.__fifos[fifo_counter].copy()
        old_score = self.__registers[register_counter].copy()

        new_score = old_score + value_in - old_max

        self.__fifos[fifo_counter]         = value_in
        self.__registers[register_counter] = new_score

        return new_score


class ScoreAccumulatorFp:
    """score accumulator, fixed-point
    """

    @property
    def q(self) -> int:
        """getter for self.__q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def p(self) -> int:
        """getter for self.__p

        Returns:
            int: value of p
        """
        return self.__p

    @property
    def nsymb(self) -> int:
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

    __fifos: np.ndarray[APUfixed]
    __registers: np.ndarray[APUfixed]

    def __init__(self, q: int, nsymb: int, bit_width: input, bit_int: int):
        self.__q         = q
        self.__p         = int(np.log2(q))
        self.__nsymb     = nsymb
        self.__mask      = q - 1
        self.__limit     = nsymb * q - 1
        self.__bit_width = bit_width
        self.__bit_int   = bit_int
        self.__fifos     = np.array([APUfixed(0, bit_width, bit_int)
                                     for _ in range(q * nsymb)])
        self.__registers = np.array([APUfixed(0, bit_width + self.p + 1, bit_int + self.p + 1)
                                     for _ in range(q)])
        self.__counter   = np.uint16(0)

    def __step_counter(self) -> np.uint16:
        curr_counter = self.__counter
        self.__counter = (curr_counter + 1) * np.uint16(curr_counter != self.__limit)
        return (curr_counter, np.bitwise_and(curr_counter, self.__mask))

    def process(self, value_in: APUfixed) -> APUfixed:
        """process a new max and update internal state

        Args:
            value_in (APUfixed): a new absolute maximum of correlation

        Returns:
            APUfixed: the resulting score
        """
        fifo_counter, register_counter = self.__step_counter()
        old_max   = self.__fifos[fifo_counter]
        old_score = self.__registers[register_counter]

        new_score = (old_score + value_in - old_max).truncate(1, False).saturate(1)

        self.__fifos[fifo_counter]         = value_in
        self.__registers[register_counter] = new_score

        return new_score

if __name__ == '__main__':
    import struct

    from matplotlib import pyplot as plt

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    GFQ    = 64
    RUNS   = 10
    NFRAME = 10
    NSYMBS = RUNS * (NFRAME + (NFRAME // 2) * 2)
    NCHIPS = NSYMBS * GFQ
    SNR    = -10.

    sigma   = np.sqrt(10**(-SNR / 10))
    sigma_c = sigma / np.sqrt(2)

    with open('data.dat', 'rb') as f:
        raw_data = f.read()
    nvalue = len(raw_data) // 4
    data   = np.array(struct.unpack(nvalue * 'f', raw_data), dtype = np.float32)

    plt.figure()
    plt.plot(data, label = 'Data (FP)')

    IN_W = 7
    IN_I = 3

    DATA_W = 13
    DATA_I = 15

    accu_flt  = ScoreAccumulator(GFQ, NFRAME)
    score_flt = np.array([accu_flt.process(value) for value in data])
    accu_fp  = ScoreAccumulatorFp(GFQ, NFRAME, DATA_W, DATA_I)
    score_fp = np.array([accu_fp.process(APUfixed(value, DATA_W, DATA_I)) for value in data])

    plt.figure()
    plt.plot(score_flt, '-', label='score_flt')
    plt.plot(score_fp, '-', label='score_fp')
    plt.show()
