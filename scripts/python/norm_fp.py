"""_summary_
"""

# pyright: basic

from multiprocessing import Pool, cpu_count

import math as m
import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt

from aptypes import APFixed, APUfixed, APComplex
from utilities import saturate, qfx


class Norm:
    @property
    def q(self):
        return self.__q

    def __init__(self, q: int, init_value: float = 2.):
        self.__q = q
        self._counter = np.uint8(0)
        self._fifo = np.ones(q, dtype=np.float32) * init_value
        self._register = np.sum(self._fifo)

    def process(self, value_in: np.complex64):
        new_value = np.real(value_in.conj() * value_in)
        counter = self._counter
        old_value = self._fifo[counter]
        self._fifo[counter] = new_value
        self._counter = np.bitwise_and(counter + 1, self.q - 1)
        old_norm = self._register
        self._register = old_norm + new_value - old_value
        return self._register


class NormFP:
    @property
    def in_width(self):
        return self._in_width

    @property
    def in_int(self):
        return self._in_int

    @property
    def in_quote(self):
        return self._in_quote

    @property
    def q(self):
        return self._q

    @property
    def p(self):
        return self._p

    @property
    def fifo_qfx(self) -> qfx:
        return self._fifo_qfx

    @property
    def reg_qfx(self) -> qfx:
        return self._reg_qfx

    @property
    def out_qfx(self) -> qfx:
        return qfx(self.__local_result.bit_width, self.__local_result.bit_int)
    
    def __init__(self, q: int, bit_width: int, bit_int: int, init_value: float | int = 1):
        self._in_width = bit_width
        self._in_int = bit_int
        self._in_quote = bit_width - bit_int
        self._q = q
        self._p = int(np.log2(q))
        self._counter = np.uint8(q - 1)
        self._fifo_qfx = qfx(bit_width, bit_int)
        self._fifo = np.array(
            [APUfixed(init_value, *self.fifo_qfx) for _ in range(q)])
        self._reg_qfx = self.fifo_qfx + self.p + 1
        self._register = APUfixed(np.sum(self._fifo), *self._reg_qfx)
        base_register_redux = (self.p // 2, int(m.ceil(self.p / 2)) + 1)
        min_quote_bits = 2
        adaptation = 0 if self._reg_qfx.bit_int - self.p // 2 + min_quote_bits < bit_width \
            else (self._reg_qfx.bit_int - self.p // 2 + min_quote_bits) - bit_width
        self._reg_rdx = (
            base_register_redux[0] + adaptation, base_register_redux[1] - adaptation)

    def __step_counter(self):
        counter = np.bitwise_and(
            self._counter + 1, self._q - 1, dtype=np.uint8)
        self._counter = counter
        return counter

    def __step_register(self, new_value):
        counter = self.__step_counter()
        old_value = self._fifo[counter]
        self._fifo[counter] = new_value
        old_norm = self._register
        # Stable
        self._register = ((old_norm + new_value) -
                          old_value).saturate(2)
        # self._register = np.sum(self._fifo).truncate(self.q)

        return self._register

    def process(self, value_in: APComplex):
        full_magn = value_in.magn().truncate(
            self.in_width - self.in_int).saturate(self.in_int + 1)
        self.__local_result = self.__step_register(full_magn).saturate(
            self._reg_rdx[0]).truncate(self._reg_rdx[1])
        return self.__local_result


if __name__ == '__main__':
    def runner(quantization: qfx, dat, gf, sigm):
        tmp_max = APFixed(0,
                          quantization.bit_width,
                          quantization.bit_int).max_value
        tmp_min = APFixed(0,
                          quantization.bit_width,
                          quantization.bit_int).min_value
        saturated_data: NDArray[np.complex64] = saturate(
            dat, tmp_max, tmp_min)  # pyright: ignore[reportAssignmentType]
        fixed_data = np.fromiter((
            APComplex(z, quantization.bit_width, quantization.bit_int).value
            for z in saturated_data),  np.complex64)
        del tmp_min, tmp_max

        norm = Norm(gf, 1. / sigm**2)
        normfp = NormFP(gf,
                        quantization.bit_width, quantization.bit_int,
                        1. / sigm**2)

        NITP = np.array([norm.process(z) for z in fixed_data])
        NITFP = np.array([normfp.process(APComplex(z, *quantization))
                         for z in fixed_data])



        print(f'[INFO] IN  >> {qfx(normfp.in_width, normfp.in_int)}')
        print(f'[INFO] OUT >> {normfp.out_qfx}')

        return [NITP, NITFP]

    pool = Pool(processes=cpu_count())

    rng = np.random.default_rng(np.random.MT19937(np.random.SeedSequence(0)))

    GFQ = 64
    RUNS = 10
    NFRAME = 10
    NSYMBS = RUNS * (NFRAME + (NFRAME // 2) * 2)
    NCHIPS = NSYMBS * GFQ
    PN = np.sign(rng.normal(size=GFQ)).astype(np.float32)
    SNR = -10.

    text = [
        'The current simulation assumes three hypotheses:',
        '  1. The channel "power" is broadly known, i.e, it is possible to define',
        '     a value `V` such that:',
        '         $\\forall x \\in X, -2^(V - 1) \\leq x \\leq 2^(V - 1)$',
        '  2. The signal is not omnipresent, i.e., the inputs are not stuck to 2^(IN_I - 1)',
        '  3. The channel is gaussian-ish, such that a value $\\sigma$ can be estimated',
        '     to compute the *mean value* of the noise.',
        '  4. The receiver gain is such that |Re(y(n))| < 1 and |Im(y(n))| < 1'
    ]
    print('\n'.join(text))

    SIGMA = np.sqrt(10**(-SNR / 10))
    SIGMA_CPX = float(SIGMA / np.sqrt(2))

    DUMMY_FRAME = np.concatenate((np.zeros(GFQ * NFRAME // 2),
                                  np.tile(PN, NFRAME),
                                  np.zeros(GFQ * NFRAME // 2)))

    data = np.array(np.sum(rng.normal(0., SIGMA_CPX, size=(NCHIPS, 2)) * (1, 1j), axis=1)
                    + np.tile(DUMMY_FRAME, RUNS),
                    dtype=np.complex64)
    IN_I = 2
    # Equivalent to channel normalization
    data_normalized = data * 2**(IN_I - 2) / max(np.max(np.abs(data.real)),
                                                 np.max(np.abs(data.imag)))
    # data_normalized = data

    norm_proc_flt = Norm(GFQ, 1. / SIGMA**2)
    norm_flt = np.array([norm_proc_flt.process(x) for x in data_normalized])

    class dummy(tuple):
        def get(self):
            return self

    result = []
    for IN_W in range(5, 16+1):
        result.append(pool.apply_async(
            runner, [qfx(IN_W, IN_I), data_normalized, GFQ, SIGMA]))
        # result.append(dummy(runner(qfx(IN_W, IN_I), data_normalized, GFQ, SIGMA)))
    figures = []
    for i, v in enumerate(result):
        norm_flt_sat, norm_fp = v.get()

        fig = plt.figure(i)
        fig.clf()
        figures.append(fig)
        norm_fp_flt = norm_fp.astype(float)

        plt.subplot(2, 1, 1)
        plt.title(f'IN_W = {i + 6}', figure=fig)
        plt.plot(norm_flt - np.mean(norm_flt), 'k:', label="Float", figure=fig)
        plt.plot(norm_flt_sat - np.mean(norm_flt_sat),
                 'b-x', label="Float Sat.", figure=fig)
        plt.plot(norm_fp_flt - np.mean(norm_fp_flt),
                 'g-x', label="FP Sat.", figure=fig)
        fig.legend()
        plt.subplot(2, 1, 2)
        plt.plot(norm_flt, 'k:', label="Float", figure=fig)
        plt.plot(norm_flt_sat,
                 'b-x', label="Float Sat.", figure=fig)
        plt.plot(norm_fp_flt,
                 'g-x', label="FP Sat.", figure=fig)
        fig.legend()

    CONDITION = True
    while CONDITION:
        FIGNO_DICT = dict(zip(range(len(figures)), range(6, 16+1)))
        print(f'Enter a figures nb ({FIGNO_DICT}), "all" or "quit"')
        value = input()
        try:
            FIGNO = int(value)
        except:
            FIGNO = -1
        if FIGNO in range(len(figures)):
            try:
                figures[FIGNO].show()
            except:
                print(f"Figure {FIGNO} has been closed already.")
        elif value == 'all':
            plt.show(block=False)
        elif value == 'quit':
            CONDITION = False
        else:
            print("Wrong input.")

    print(__file__ + ': ok')
