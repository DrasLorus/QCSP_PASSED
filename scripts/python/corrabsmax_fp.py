"""TSCorrAbsMax Fixed-point implementation
"""

# pyright: basic

from multiprocessing import Pool, cpu_count 

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt 

from aptypes import APComplex, APUfixed, APFixed

from utilities import extract_data, saturate
from correlation_fp import TimeSlidingCorrelatorFP

from corrabsmax import TSCorrAbsMax

class TSCorrAbsMaxFP:
    """extract the absolute max of correlation, fixed-point version
    """

    @property
    def p(self):
        return self.__correlator.p

    @property
    def bit_width(self):
        return self.__correlator.bit_width

    @property
    def bit_int(self):
        return self.__correlator.bit_int

    def __init__(self, gf_q: int, pn: NDArray[np.float32], bit_width: int = 16, bit_int: int = 6) -> None:
        self.__correlator = TimeSlidingCorrelatorFP(
            gf_q, pn, bit_width, bit_int)

    @classmethod
    def process(cls, *_):
        raise RuntimeError(
            'Fixed-point version only support square output. Use process_sqr.')

    def process_sqr(self, value_in: APComplex):
        """process value_in through the complete correlation process

        Args:
            value_in (APComplex): values to process

        Returns:
            np.floating: resulting absolute maxima of correlation.
        """
        correlation = self.__correlator.process(value_in)
        full_max_values: list[APUfixed] = [z.magn() for z in correlation]
        full_max_value = max(full_max_values)
        #full_max_values: APUfixed = max(
        #    tuple(map(lambda z: z.magn(), self.__correlator.process(value_in))))
        return full_max_value.saturate(int(self.p) + 1).truncate(self.bit_width + 1)


def play(var_names: list[str], qfx = [16, 6]):
    GFQ = 64
    NFR = 60
    COUNT = 5
    RUNS = 1
    NBR = GFQ * NFR * COUNT * RUNS

    KEY = var_names[2]

    pn = extract_data('data/parameters_20210903.mat',
                      ['PN64'])['PN64'].reshape((64))
    dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names[:2])

    DATA: NDArray[np.complex64] = dict_data[var_names[0]
                                            ][0][:NBR].astype(np.complex64)
    CABSMAX_TRUE = dict_data[var_names[1]][0][:NBR].astype(np.float32)
    PN = pn.astype(np.float32)

    IN_W = qfx[0]
    IN_I = qfx[1]

    TMP_MAX = APFixed(0, IN_W, IN_I).max_value
    TMP_MIN = APFixed(0, IN_W, IN_I).min_value
    SATURATED_DATA: NDArray[np.complex64] = saturate(
        DATA, TMP_MAX, TMP_MIN)  # pyright: ignore[reportAssignmentType]
    del TMP_MIN, TMP_MAX

    FIXED_DATA = np.fromiter((APComplex(z, IN_W, IN_I)
                              for z in SATURATED_DATA), dtype=APUfixed)

    cabsmaxfl = TSCorrAbsMax(GFQ, PN)
    cabsmaxfl_fp = np.array([cabsmaxfl.process_sqr(z.value)
                            for z in FIXED_DATA])

    cabsmaxfp = TSCorrAbsMaxFP(GFQ, PN, IN_W, IN_I)
    cabsmaxfp_fp = np.array([cabsmaxfp.process_sqr(z).value for z in FIXED_DATA])

    print(f'-- Process "{KEY} [{IN_W}]" ok.')

    return (CABSMAX_TRUE, cabsmaxfl_fp, cabsmaxfp_fp, KEY, IN_W)

def main():
    var_names: dict[str, list[str]] = dict()
    var_names['sqr_raw_m10dB'] = ['data_input_m10dB_w1_q64_N60_0_n30',
                                  'cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30']
    var_names['sqr_raw_infdB'] = ['data_input_infdB_w1_q64_N60_0_n10',
                                  'cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10']
    GFQ = 64

    IN_I = 6
    IN_Ws = range(8, 17)
    pool = Pool(min(cpu_count(), len(IN_Ws)))

    results = []
    for IN_W in IN_Ws:
        #for key,var in var_names.items():
        #    results.append(pool.apply_async(play, [var + [key], [IN_W, IN_I]]))
        key = 'sqr_raw_m10dB'
        results.append(pool.apply_async(play, [var_names[key] + [key], [IN_W, IN_I]]))
    for i,v in enumerate(results):
        values = v.get()
        KEY = values[3]
        IN_W = values[4]
        plt.figure()
        plt.title(f'{KEY}: from chip {10*GFQ} to {16*GFQ}')
        plt.plot(values[0][10*GFQ:16*GFQ], '-o', label='True Float')
        plt.plot(values[1][10*GFQ:16*GFQ], '-d', label='Float on Fixed')
        plt.plot(values[2][10*GFQ:16*GFQ], '-x', label='True Fixed')
        plt.legend()

        plt.figure()
        plt.title(f'{KEY}: complete')
        plt.plot(values[0][GFQ:], '-o', label='True Float')
        plt.plot(values[1][GFQ:], '-d', label='Float on Fixed')
        plt.plot(values[2][GFQ:], '-x', label='True Fixed')
        plt.legend()

        print(f'-- Figure {i+1}/{len(results)}  ({KEY} [{IN_W}]) ok.')
            
    plt.show(block=False)
    #_ = [play(variables + [key], qfx) for key, variables in var_names.items()]


if __name__ == '__main__':
    main()
    input('Hit any key to end...')

    print(__file__ + ': ok')
