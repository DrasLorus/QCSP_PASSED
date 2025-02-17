import numpy as np

class QSpannedSequentialAdder:
    """Iterative filter involved in Time sliding-based correlations
    """
    @property
    def q(self) -> int:
        """value of q

        Returns:
            int: value of q
        """
        return self.__q

    @property
    def counter(self) -> np.uint16:
        """getter for the counter

        Returns:
            np.uint16: value of the counter
        """
        return self.__counter

    def __init__(self, q: int):
        self.__q = np.uint16(q)
        self.__counter = np.uint16(0)
        self.__fifo = np.zeros(q, dtype=np.complex64)

    def __step_counter(self) -> np.uint16:
        """increment and return the value of the counter

        Returns:
            np.uint16: value of the stepped counter
        """
        counter = self.__counter
        self.__counter = np.bitwise_and(counter + 1, self.__q - 1).astype(np.uint16)
        return counter

    def process(self, value_in : np.complex64) -> np.complex64:
        """process value_in through the iterative filter

        Args:
            value_in (np.complex64): new value

        Returns:
            np.complex64: resulting value
        """
        counter = self.__step_counter()
        old_value = self.__fifo[counter]
        self.__fifo[counter] = value_in
        return value_in - old_value
