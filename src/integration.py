import numpy as np


def point_matching_integral(a, x, y):
    zpa = complex(x, y) + a
    zma = complex(x, y) - a
    term_1 = zma * (np.log(zma) - 1)
    term_2 = zpa * (np.log(zpa) - 1)
    total = term_1 - term_2
    return total.real


def galerkin_integral(a, x1, y1, x2, y2):
    z1 = complex(x1, y1)
    z2 = complex(x2, y2)
    z1ma = complex(x1, y1) - a
    z2ma = complex(x2, y2) - a
    z1pa = complex(x1, y1) + a
    z2pa = complex(x2, y2) + a

    def func(z):
        if z == 0:
            return 0
        else:
            return 0.5 * z ** 2 * (np.log(z) - 1.5)

    p1 = np.array((x1, y1))
    p2 = np.array((x2, y2))
    d = p2-p1
    norm_d = np.linalg.norm(d)
    if norm_d == 0:
        raise ValueError('points cannot be identical')
    else:
        d_hat = d / norm_d
    l_hat = complex(d_hat[0], d_hat[1])
    #l_hat = (z2 - z1) / abs(z2 - z1)

    term1 = func(z2ma)
    term2 = func(z2pa)
    term3 = func(z1ma)
    term4 = func(z1pa)

    k = l_hat.conjugate() * ((term1 - term2) - (term3 - term4))

    return k.real


