"""utility functions
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.random as rd


import h5py as h5

def saturate(z: complex, max_value: float, min_value: float) -> complex:
    """saturate real and imag part to min and max values

    Args:
        z (complex): the initial complex
        max_value (float): maximum real value
        min_value (float): minimum real value

    Returns:
        complex: saturated value
    """
    return min(max(z.real, min_value), max_value) + 1j*min(max(z.imag, min_value), max_value)

@dataclass(frozen=True)
class qfx:
    bit_width: int
    bit_int: int

    @property
    def bit_quote(self) -> int:
        return self.bit_width - self.bit_int

    def __add__(self, value: list[int] | tuple[int] | qfx) -> qfx:
        if isinstance(value, qfx):
            return qfx(self.bit_width + value.bit_width, self.bit_int + value.bit_int)
        if isinstance(value, (int, np.integer)):
            return qfx(self.bit_width + value, self.bit_int + value)
        return qfx(self.bit_width + value[0], self.bit_int + value[1])
    def __iter__(self):
        return iter((self.bit_width, self.bit_int))
    def __str__(self):
        return f'QFX :: Width = {self.bit_width} | Int. = {self.bit_int} | Quote = {self.bit_quote}'

def generate_data(gf_q: int,
                  runs: int,
                  snr: float,
                  pn: np.ndarray[np.float32] = None,
                  rng: int | rd.Generator = None) -> np.ndarray[np.complex64]:
    """_summary_

    Args:
        gf_q (int): _description_
        nbr (int): _description_
        pn (numpy.ndarray[numpy.float32]): _description_
        snr (float): _description_
        rng (int | numpy.random.Generator, optional): _description_. Defaults to None.

    Returns:
        numpy.ndarray[numpy.complex64]: _description_
    """
    if rng is None:
        rng = rd.default_rng(rd.MT19937(rd.SeedSequence(rd.random_integers(0, 2**32))))
    elif isinstance(rng, int):
        rng = rd.default_rng(rd.MT19937(rd.SeedSequence(rng)))
    elif isinstance(rng, rd.Generator):
        pass
    else:
        raise TypeError(f"rng cannot be a {type(rng)}")

    spl_nb = runs * gf_q
    pn_seq = np.sign(rd.randn(gf_q)).astype(np.float32) if pn is None else pn

    sigma   = np.sqrt(10**(-snr / 10))
    sigma_c = sigma / np.sqrt(2)


    data = np.array(np.sum(rng.normal(0.,sigma_c, size=(spl_nb, 2)) * (1, 1j), axis=1)
                 + np.tile(pn_seq, spl_nb // gf_q),
                dtype=np.complex64)
    return data

def extract_variable(h5_file: h5.File, var_name: str) -> np.ndarray[np.float32 | np.complex64]:
    """_summary_

    Args:
        h5_file (h5.File): _description_
        var_name (str): _description_

    Returns:
        np.ndarray[np.float32 | np.complex64]: _description_
    """
    _local = h5_file[var_name]
    _dtype = _local[0].dtype
    if _dtype == 'float32':
        local_data = np.array(_local)
    elif _dtype == [('real', '<f4'), ('imag', '<f4')]:
        _re = _local[:]['real']
        _im = _local[:]['imag']
        local_data = _re + 1j * _im
    else:
        raise RuntimeError(f'Variable "{var_name}" has invalide type "{_dtype}".')
    return local_data

def extract_data(filepath: Path | str, var_names: list[str]) -> dict[str:np.ndarray]:
    """_summary_

    Args:
        file (Path | str): _description_

    Returns:
        numpy.ndarray: _description_
    """
    file = h5.File(filepath)

    extracted = dict()
    for var in var_names:
        extracted[var] = extract_variable(file, var)

    file.close()
    del(file)
    return extracted
