"""utility functions
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np

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
