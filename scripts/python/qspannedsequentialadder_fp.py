from aptypes import APComplex

import numpy as np

class QSpannedSequentialAdderFP:
    """perform the 'iterative accumulation'

    """
    @property
    def q(self) -> np.uint16:
        """gf cardinal

        Returns:
            np.uint16: value of q
        """
        return self.__q

    @property
    def counter(self) -> np.uint16:
        """internal counter

        Returns:
            np.uint16: the value
        """
        return self.__counter

    @property
    def bit_width(self) -> int:
        """input quantization

        Returns:
            int: the value
        """
        return self._in_width

    @property
    def bit_int(self):
        """input quantization, int-part

        Returns:
            int: the value
        """
        return self._in_int

    @property
    def bit_quote(self):
        """input quantization, quote-part

        Returns:
            int: the value
        """
        return self._in_quote

    def __init__(self, q: int, bit_width = 16, bit_int = 4):
        self._in_width = bit_width
        self._in_int   = bit_int
        self._in_quote = bit_width - bit_int
        self.__q       = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo    = np.zeros((q, 2), dtype=int)

    def __step_counter(self) -> np.uint16:
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1)
        return counter

    def __update_fifo(self, value_in: APComplex) -> APComplex:
        """update the fifo with the new value and get the oldest one

        Args:
            value_in (APComplex): new value

        Returns:
            APComplex: oldest value in the fifo
        """
        counter              = self.__step_counter()
        old_value            = self.__fifo[counter].copy()
        self.__fifo[counter] = (value_in.real.raw, value_in.imag.raw)
        return APComplex(old_value, self.bit_width, self.bit_int)

    def process(self, value_in: APComplex) -> APComplex:
        """process value_in through the acumulative filter

        Args:
            value_in (APComplex): new value in the filter

        Returns:
            APComplex: the newest value minus the oldest one
        """
        old_value = self.__update_fifo(value_in)
        return value_in - old_value
