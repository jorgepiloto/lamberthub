import numpy as np


def _get_alpha0(a, s):
    if a > 0:
        alpha0 = 2 * np.arcsin(np.sqrt((s) / (2 * a)))
    else:
        alpha0 = 2 * np.arcsinh(np.sqrt((s) / (-2 * a)))
    return alpha0


def _get_beta0(a, c, s):
    if a > 0:
        beta0 = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
    else:
        beta0 = 2 * np.arcsinh(np.sqrt((s - c) / (-2 * a)))
    return beta0


def _get_orbit_parameter(a, r1_norm, r2_norm, c, s, alpha, beta):

    if a >= 0:
        p = (
            (4 * a * (s - r1_norm) * (s - r2_norm))
            * np.sin((alpha + beta) / 2) ** 2
            / (c ** 2)
        )
    else:
        p = (
            (-4 * a * (s - r1_norm) * (s - r2_norm))
            * np.sinh((alpha + beta) / 2) ** 2
            / (c ** 2)
        )
    return p
