"""Comparison of correlation method approach
"""
from argparse import ArgumentParser

import numpy as np
from matplotlib import pyplot as plt

from utilities import generate_data,extract_data
from qspannedsequentialadder import QSpannedSequentialAdder

class TimeSlidingCorrelator(QSpannedSequentialAdder):
    """Time sliding-based correlator, floating-point implementation 
    """
    @property
    def pn(self) -> np.ndarray[np.float32]:
        """the PN sequence

        Returns:
            np.ndarray[np.float32]: the PN sequence
        """
        return self.__pn

    def __rotate_pn(self) -> None:
        """rotate the PN sequence in place
        """
        self.__pn = np.roll(self.__pn, -1)

    def __init__(self, q: int, pn: np.ndarray):
        super().__init__(q)
        self.__pn = pn
        self.__registers = np.zeros(q, dtype=np.complex64)

    def __permute(self, correlations: np.ndarray[np.complex64]):
        """perform the permutation

        Args:
            correlations (np.ndarray[np.complex64]): correlation array

        Returns:
            np.ndarray[np.complex64]: permuted correlation array
        """
        cnt = self.counter
        return np.flip(np.roll(correlations, int(cnt) - 1))


    def __pn_correlate(self, value_in: np.complex64):
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        new_values       = value_in * self.pn
        old_correlations = self.__registers
        new_correlations = old_correlations + new_values
        self.__registers = new_correlations
        self.__rotate_pn()
        return new_correlations

    def process(self, value_in: np.complex64):
        """correlate a new stepped value with the PN sequence

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        return self.__pn_correlate(super().process(value_in))

    def process_permuted(self, value_in : np.complex64) -> np.ndarray[np.complex64]:
        """retrieve the permuted correlation, matching FFT-based ones

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: permuted correlation vector
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
    def pn(self) -> np.ndarray[np.float32]:
        """getter for the PN sequence

        Returns:
            np.ndarray[np.float32]: the sequence
        """
        return np.fft.ifft(self.__cfpn.conj()).real.astype(np.float32)

    def __init__(self, q: int, pn: np.ndarray[np.float32]):
        self.__q = np.uint16(q)
        self.__cfpn = np.fft.fft(pn).astype(np.complex64).conj()
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def process(self, value_in : np.complex64) -> np.ndarray[np.complex64]:
        """perform the correlation using FFT method

        Args:
            value_in (np.complex64): new value

        Returns:
            np.ndarray[np.complex64]: correlation vector
        """
        self.__fifo = np.concatenate((self.__fifo[1:], [value_in]))
        return np.fft.ifft(np.fft.fft(self.__fifo) * self.__cfpn)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--data', type=str, default='generate')
    args = parser.parse_args()

    GFQ  = 64
    RUNS = 10**3
    PN   = np.sign(np.random.randn(GFQ)).astype(np.float32)
    SNR  = -10
    SEED = 0

    if args.data == "generate":
        data = generate_data(GFQ, RUNS, SNR, PN, SEED)
    elif args.data == "extract":
        var_names = ['data_input_m10dB_w1_q64_N60_0_n30', 'cabs_max_raw_m10dB_w1_q64_N60_0_n30']
        dict_data = extract_data('data/test_data_w1_nofreq.mat', var_names)
        data = dict_data['data_input_m10dB_w1_q64_N60_0_n30'][0][:RUNS*GFQ]

    fft_corr_proc = FftCorrelator(GFQ, PN)
    ts_corr_proc  = TimeSlidingCorrelator(GFQ, PN)

    fft_corr = np.array([fft_corr_proc.process(z) for z in data])
    ts_corr  = np.array([ts_corr_proc.process_permuted(z) for z in data])

    total_size = len(ts_corr)

    plt.figure()
    plt.plot(fft_corr[10*GFQ:16*GFQ, 0].real, label = 'FFT')
    plt.plot(ts_corr[10*GFQ:16*GFQ, 0].real, label = 'TS')


    plt.figure()
    plt.plot(fft_corr[12*GFQ, :].real, label = 'FFT')
    plt.plot(ts_corr[12*GFQ, :].real, label = 'TS')

    plt.figure()
    plt.plot([np.max(np.abs(fft_corr[i-GFQ:i])) for i in range(GFQ, total_size)], label = 'FFT')
    plt.plot([np.max(np.abs(ts_corr[i-GFQ:i])) for i in range(GFQ, total_size)], label = 'TS')

    if args.data == "extract":
        plt.figure()
        plt.plot(dict_data['cabs_max_raw_m10dB_w1_q64_N60_0_n30'][0][GFQ:GFQ:(RUNS+1)*GFQ], label = 'TS True')
        plt.plot([np.max(np.abs(ts_corr[i-GFQ:i])) for i in range(GFQ, total_size)], label = 'TS')

    plt.show()

    print(__file__ + ': ok')
