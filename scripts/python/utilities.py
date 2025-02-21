"""utility functions
"""

# pyright: basic, reportAssignmentType=false

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy import random as rd
from numpy.typing import NDArray

import h5py as h5


def saturate(z: complex | NDArray[np.complexfloating],
             max_value: float,
             min_value: float) -> complex | NDArray[np.complexfloating]:
    """saturate real and imag part to min and max values

    Args:
        z (complex): the initial complex
        max_value (float): maximum real value
        min_value (float): minimum real value

    Returns:
        complex: saturated value
    """
    return np.clip(z.real, min_value, max_value) + 1j * np.clip(z.imag, min_value, max_value)


@dataclass(frozen=True)
class qfx:
    bit_width: int
    bit_int: int

    @property
    def bit_quote(self) -> int:
        return self.bit_width - self.bit_int

    def __add__(self, value: int | tuple[int, int] | qfx) -> qfx:
        if isinstance(value, qfx):
            return qfx(self.bit_width + value.bit_width, self.bit_int + value.bit_int)
        if isinstance(value, (int, np.integer)):
            return qfx(self.bit_width + value, self.bit_int + value)
        return qfx(self.bit_width + value[0], self.bit_int + value[1])

    def __iter__(self):
        return iter((self.bit_width, self.bit_int))

    def __str__(self) -> str:
        return f'QFX :: Width = {self.bit_width} | Int. = {self.bit_int} | Quote = {self.bit_quote}'


def generate_data(gf_q: int,
                  nbr: int,
                  snr: float,
                  pn: NDArray[np.float32] | None = None,
                  rng: int | rd.Generator | None = None) -> NDArray[np.complex64]:
    """ generate random values mixing pn and gaussian noise

    Args:
        gf_q (int): Galois Field size.
        nbr (int): Number of "symbols" of gf_q samples to generate.
        pn (numpy.ndarray[numpy.float32]): A PN sequence or `None` to auto generate it. Defaults to None.
        snr (float): Targeted signal-to-noise ratio. 
        rng (int | numpy.random.Generator, optional): A random number generator or a seed. Defaults to a randomly seeded default_rng based on an MT19937.

    Returns:
        numpy.ndarray[numpy.complex64]: generated data
    """
    if rng is None:
        rng = rd.default_rng(rd.MT19937(
            rd.SeedSequence(rd.random_integers(0, 2**32))))
    elif isinstance(rng, int):
        rng = rd.default_rng(rd.MT19937(rd.SeedSequence(rng)))

    spl_nb = nbr * gf_q
    pn_seq = np.sign(rd.randn(gf_q)).astype(np.float32) if pn is None else pn

    SIGMA: float = np.sqrt(10**(-snr / 10))
    SIGMA_C: float = SIGMA / np.sqrt(2)

    data = np.array(np.sum(rng.normal(0., SIGMA_C, size=(spl_nb, 2)) * (1, 1j), axis=1)
                    + np.tile(pn_seq, spl_nb // gf_q),
                    dtype=np.complex64)
    return data


def extract_variable(h5_file: h5.File, var_name: str) -> NDArray[np.float32 | np.complex64]:
    """ extract a variable from an opened HDF5 or MAT-7.3 file 

    Args:
        h5_file (h5py.File): h5py file descripto
        var_name (str): variable name in the file

    Returns:
        NDArray[np.float32 | np.complex64]: extracted variable
    """

    _local: h5.Dataset = h5_file[var_name]
    _dtype = _local[0].dtype
    if _dtype == 'float32' or _dtype == 'float64':
        local_data = np.array(_local)
    elif _dtype == [('real', '<f4'), ('imag', '<f4')]:
        _re: NDArray[np.float32] = _local[:]['real']
        _im: NDArray[np.float32] = _local[:]['imag']
        local_data = (_re + 1j * _im).astype(np.complex64)
    else:
        raise RuntimeError(
            f'Variable "{var_name}" has invalide type "{_dtype}".')
    return local_data


def extract_data(filepath: Path | str, var_names: list[str]) -> dict[str, NDArray[np.float32 | np.complex64]]:
    """ handle extraction of several variables from an HDF5 or MAT-7.3 file

    Args:
        file (Path | str): path to the file
        var_names (list[str]): list of variables to extract

    Returns:
        NDArray: array of variables
    """
    file = h5.File(filepath, 'r')

    extracted: dict[str, NDArray[np.float32 | np.complex64]] = dict()
    for var in var_names:
        extracted[var] = extract_variable(file, var)

    file.close()
    del file
    return extracted
