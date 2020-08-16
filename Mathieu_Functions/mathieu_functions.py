"""
defines a class where the Mathieu functions ce and se are defined.
"""
import numpy as _np
from eig_system import matrix_system, eig_pairs


class mathieu_functions:

    def __init__(self, q, N, type, period):
        self._q = q
        self._N = N
        self._type = type
        self._period = period

    @classmethod
    def ce_even(
        cls,
        q,
        x,
        N,
        type='even',
        period='one',
        debug=False
    ):
        """Cosine elliptic function ce_{2n}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `a`, and the domain.
        """
        a, As = eig_pairs(matrix_system(q, N))
        vals = {}
        for n in range(N):
            vals.update({'a' + str(2 * n): a[n]})
            if debug:
                vals.update({'A' + str(2 * n): As[n, :]})
            ce = As[n, 0]
            for k in range(1, N):
                ce = ce + As[n, k] * _np.cos((2 * k) * x)
            vals.update({'ce' + str(2 * n): ce})
        return vals
