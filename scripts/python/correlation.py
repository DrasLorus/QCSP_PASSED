"""Comparison of correlation method approach
"""

# pyright: basic

from argparse import ArgumentParser

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt

from utilities import generate_data, extract_data
from qspannedsequentialadder import QSpannedSequentialAdder

class TimeSlidingCorrelator:
    """Time sliding-based correlator, floating-point implementation
    """

    @property
    def pn(self) -> NDArray[np.float32]:
        """the PN sequence

        Returns:
            NDArray[np.float32]: the PN sequence
        """
        return self.__pn

    def __rotate_pn(self) -> None:
        """rotate the PN sequence in place
        """
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: NDArray[np.float32]):
        self.__adder = QSpannedSequentialAdder(q)
        self.__root_pn = pn.copy()
        self.__pn = pn.copy()
        self.__registers = np.zeros(q, dtype=np.complex64)

    def reset(self) -> None:
        self.__adder.reset()
        self.__pn = self.__root_pn.copy()
        self.__registers = np.zeros(self.__adder.q, dtype=np.complex64)

    def __permute(self, correlations: NDArray[np.complex64]):
        """perform the permutation

        Args:
            correlations (NDArray[np.complex64]): correlation array

        Returns:
            NDArray[np.complex64]: permuted correlation array
        """
        cnt = self.__adder.counter
        return np.flip(np.roll(correlations, int(cnt) - 1))


    def __pn_correlate(self, value_in: np.complex64):
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            NDArray[np.complex64]: correlation vector
        """
        new_values: NDArray[np.complex64]       = value_in * self.pn
        old_correlations: NDArray[np.complex64] = self.__registers
        new_correlations: NDArray[np.complex64] = old_correlations + new_values
        self.__registers: NDArray[np.complex64] = new_correlations
        self.__rotate_pn()
        return new_correlations

    def process(self, value_in: np.complex64) -> NDArray[np.complex64]:
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            NDArray[np.complex64]: correlation vector
        """
        return self.__pn_correlate(self.__adder.process(value_in))

    def process_permuted(self, value_in : np.complex64) -> NDArray[np.complex64]:
        """retrieve the permuted correlation, matching FFT-based ones

        Args:
            value_in (np.complex64): new value

        Returns:
            NDArray[np.complex64]: permuted correlation vector
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
    def pn(self) -> NDArray[np.float32]:
        """getter for the PN sequence

        Returns:
            NDArray[np.float32]: the sequence
        """
        return np.fft.ifft(self.__cfpn.conj()).real.astype(np.float32)

    def __init__(self, q: int, pn: NDArray[np.float32]):
        self.__q = np.uint16(q)
        self.__cfpn = np.fft.fft(pn).astype(np.complex64).conj()
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def process(self, value_in : np.complex64) -> NDArray[np.complex64]:
        """perform the correlation using FFT method

        Args:
            value_in (np.complex64): new value

        Returns:
            NDArray[np.complex64]: correlation vector
        """
        self.__fifo = np.concatenate((self.__fifo[1:], [value_in]))
        return np.fft.ifft(np.fft.fft(self.__fifo) * self.__cfpn).astype(np.complex64)

def main():
    parser = ArgumentParser()
    parser.add_argument('-d', '--data', type=str,
        default='generate', choices=['generate', 'extract'])
    args = parser.parse_args()

    GFQ  = 64
    RUNS = 10**3
    SNR  = -10
    SEED = 0

    DATA_METHOD: str = args.data 

    EXTRACTED = DATA_METHOD == 'extract'
    GENERATED = DATA_METHOD == 'generate'

    if GENERATED:
        pn = np.sign(np.random.randn(GFQ))
        data = generate_data(GFQ, RUNS, SNR, pn, SEED)
    elif EXTRACTED:
        pn = extract_data('data/parameters_20210903.mat', ['PN64'])['PN64'].reshape((64))
        var_names = ['data_input_m10dB_w1_q64_N60_0_n30', 'cabs_max_raw_m10dB_w1_q64_N60_0_n30']
        dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
        data = dict_data['data_input_m10dB_w1_q64_N60_0_n30'][0][:RUNS*GFQ]
    DATA = data.astype(np.complex64)
    PN = pn.astype(np.float32)

    fft_corr_proc = FftCorrelator(GFQ, PN)
    ts_corr_proc  = TimeSlidingCorrelator(GFQ, PN)

    fft_corr = np.array([fft_corr_proc.process(z) for z in DATA])  
    ts_corr  = np.array([ts_corr_proc.process_permuted(z) for z in DATA])  

    total_size = len(ts_corr)

    plt.figure()
    plt.plot(fft_corr[10*GFQ:16*GFQ, 0].real, label = 'FFT')
    plt.plot(ts_corr[10*GFQ:16*GFQ, 0].real, label = 'TS')
    plt.legend()

    plt.figure()
    plt.plot(fft_corr[12*GFQ, :].real, label = 'FFT')
    plt.plot(ts_corr[12*GFQ, :].real, label = 'TS')
    plt.legend()

    plt.figure()
    plt.plot(np.max(np.abs(fft_corr), axis=1), label = 'FFT')
    plt.plot(np.max(np.abs(ts_corr), axis=1),  label = 'TS')
    if EXTRACTED:
        plt.plot(dict_data[var_names[1]][0][:total_size],
            label = 'TS True')
    plt.legend()

    plt.show()

    print(__file__ + ': ok')


if __name__ == '__main__':
    main()