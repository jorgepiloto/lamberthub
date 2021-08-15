import numpy as np


def _get_alpha0(a, s):
    alpha0 = 2 * np.arcsin(np.sqrt((s)/(2 * a)))
    return alpha0

def _get_beta0(a, c, s):
    beta0 = 2 * np.arcsin(np.sqrt((s - c)/(2 * a)))
    return beta0
