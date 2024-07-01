"""utility functions
"""

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
