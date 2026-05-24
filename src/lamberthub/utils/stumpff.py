"""Module containing the implementation of Stumpff functions.

These set of functions usually appear in a plethora of equations within the
subject of astrodynamics.
"""

from math import gamma

from numba import njit as jit
import numpy as np


@jit(cache=True)
def c2(psi):
    r"""Second Stumpff function.

    For positive arguments:

    .. math::

        c_2(\psi) = \frac{1 - \cos{\sqrt{\psi}}}{\psi}

    """
    eps = 1.0
    if psi > eps:
        res = (1 - np.cos(np.sqrt(psi))) / psi
    elif psi < -eps:
        res = (np.cosh(np.sqrt(-psi)) - 1) / (-psi)
    else:
        res = 1.0 / 2.0
        delta = (-psi) / gamma(2 + 2 + 1)
        k = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / gamma(2 * k + 2 + 1)

    return res


@jit(cache=True)
def c3(psi):
    r"""Third Stumpff function.

    For positive arguments:

    .. math::

        c_3(\psi) = \frac{\sqrt{\psi} - \sin{\sqrt{\psi}}}{\sqrt{\psi^3}}

    """
    eps = 1.0
    if psi > eps:
        sqrt_psi = np.sqrt(psi)
        res = (sqrt_psi - np.sin(sqrt_psi)) / (psi * sqrt_psi)
    elif psi < -eps:
        sqrt_neg_psi = np.sqrt(-psi)
        res = (np.sinh(sqrt_neg_psi) - sqrt_neg_psi) / (-psi * sqrt_neg_psi)
    else:
        res = 1.0 / 6.0
        delta = (-psi) / gamma(2 + 3 + 1)
        k = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / gamma(2 * k + 3 + 1)

    return res
