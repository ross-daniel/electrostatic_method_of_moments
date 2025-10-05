import numpy as np
import constants


def c_prime(w, h, eps_reff):
    v_p = constants.c0 / np.sqrt(eps_reff)

    wh_ratio = w / h

    if wh_ratio > 1:
        Z0 = (constants.eta_0 / np.sqrt(eps_reff)) * ((wh_ratio + 1.393 + 0.667 * np.log(wh_ratio + 1.444)) ** -1)
    else:
        Z0 = (constants.eta_0 / (2 * np.pi * np.sqrt(eps_reff))) * np.log((8 * h / w) + w / (4 * h))

    C_prime = 1 / (Z0 * v_p)
    return C_prime
