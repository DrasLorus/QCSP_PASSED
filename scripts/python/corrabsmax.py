"""TSCorrAbsMax floating-point implementation
"""

# pyright: basic

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt

from utilities import extract_data
from correlation import TimeSlidingCorrelator


class TSCorrAbsMax:
    """extract the absolute max of correlation
    """

    def __init__(self, gf_q: int, pn: NDArray[np.float32]):
        self.__correlator = TimeSlidingCorrelator(gf_q, pn)

    def process_sqr(self, value_in: np.complex64) -> np.float32:
        """process value_in through the complete correlation process

        Args:
            value_in (APComplex): new value to process

        Returns:
            np.floating: resulting absolute maxima of correlation.
        """
        return np.max(np.abs(self.__correlator.process(value_in))**2)

    def process(self, value_in: np.complex64) -> np.float32:
        """process value_in through the complete correlation process

        Args:
            value_in (APComplex): new value to process

        Returns:
            np.floating: resulting absolute maxima of correlation.
        """
        return np.max(np.abs(self.__correlator.process(value_in)))


def play(var_names: list[str]):
    GFQ = 64
    NFR = 60
    COUNT = 5
    RUNS = 3
    NBR = GFQ * NFR * COUNT * RUNS

    KEY = var_names[2]

    pn = extract_data('data/parameters_20210903.mat',
                      ['PN64'])['PN64'].reshape((64))
    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names[:2])

    DATA = dict_data[var_names[0]][0][:NBR].astype(np.complex64)
    CABSMAX_TRUE = dict_data[var_names[1]][0][:NBR].astype(np.float32)
    PN = pn.astype(np.float32)

    cabsmax = TSCorrAbsMax(GFQ, PN)
    if '_sqr_' in var_names[1]:
        cabsmax_result = np.array([cabsmax.process_sqr(z) for z in DATA])
    else:
        cabsmax_result = np.array([cabsmax.process(z) for z in DATA])

    plt.figure()
    plt.title(f'{KEY}: from chip {10*GFQ} to {16*GFQ}')
    plt.plot(cabsmax_result[10*GFQ:16*GFQ], label='Float')
    plt.plot(CABSMAX_TRUE[10*GFQ:16*GFQ], label='True')
    plt.legend()

    plt.figure()
    plt.title(f'{KEY}: complete')
    plt.plot(cabsmax_result, label='Float')
    plt.plot(CABSMAX_TRUE,  label='True')
    plt.legend()

    plt.show(block=False)


if __name__ == '__main__':
    var_names: dict[str, list[str]] = dict()
    var_names['m10dB'] = ['data_input_m10dB_w1_q64_N60_0_n30',
                          'cabs_max_raw_m10dB_w1_q64_N60_0_n30']
    var_names['sqr_m10dB'] = ['data_input_m10dB_w1_q64_N60_0_n30',
                              'cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30']
    var_names['infdB'] = ['data_input_infdB_w1_q64_N60_0_n10',
                          'cabs_max_raw_infdB_w1_q64_N60_0_n10']
    var_names['sqr_infdB'] = ['data_input_infdB_w1_q64_N60_0_n10',
                              'cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10']
    _ = [play(var + [key]) for key, var in var_names.items()]
    input('Hit any key to end...')

    print(__file__ + ': ok')
