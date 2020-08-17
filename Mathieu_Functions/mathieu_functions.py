"""
defines a class where the Mathieu functions ce and se are defined.
"""
from eig_system import matrix_system, eig_pairs
from math import factorial
import numpy as _np


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
        vals = {}
        for n in range(N):
            a, As = eig_pairs(matrix_system(q[0], N))
            a = [a[n]]  # makes a list of the nth eigenvalue
            As = As[_np.newaxis, :, n]
            for k in range(1, len(q)):
                an, nAs = eig_pairs(matrix_system(q[k], N))
                a.append(an[n])
                nAs = nAs[_np.newaxis, :, n]
                As = _np.append(As, nAs, axis=0)
            As = Fcoeffs(As, n)
            vals.update({'a' + str(2 * n): _np.array(a)})
            vals.update({'A' + str(2 * n): As})
        # initialize the coefficients (q=0)
        for n in range(N // 4 + 1):
            vals['ce' + str(2 * n)] = 0
            terms = [_np.cos((2 * k) * x) * (vals['A' + str(2 * n)][0, k])
                     for k in range(N // 4 + 1)]
            vals.update({'ce' + str(2 * n): (vals['ce' + str(2 * n)]) +
                         + _np.sum(terms, axis=0)})
            vals.update({'ce' + str(2 * n):
                         vals['ce' + str(2 * n)][_np.newaxis, :]})
        for i in range(1, len(q)):
            for n in range(N // 4 + 1):
                terms = [_np.cos((2*k)*x)*(vals['A'+str(2*n)][i-1, k])
                         for k in range(N // 4 + 1)]
                ce = vals['ce' + str(2 * n)][i - 1, :] + _np.sum(terms, axis=0)
                ce = ce[_np.newaxis, :]
                ce = _np.append(vals['ce' + str(2 * n)], ce, axis=0)
                vals.update({'ce' + str(2 * n): ce})
        return vals


def Fcoeffs(As, n=0, q=0.001):
    """ Returns the Fourier coefficient of the Mathieu functions for given
    parameter q. Makes sure the coefficients are continuous (in q). Numerical
    routines for estimating eigenvectors might converge in different signs
    for the eigenvectors for different (neighboring) values of q.
        Input:
            As: 1d array. Eigenvector shape(As)=Nq, N, as a function of q and
                containing N entries, each associated with a Fourier
                coefficient.
            n: int. Eigenvalue index. If n=0, eigenvalue is a0. n=1, eigenvalue
            is a2.
            q: float, real or imag. Default is q=0.01, real. Must be small.
        Output:
            Corrected Eigenvector with same shape as original
    """
    # Estimate limiting value for small q (pos or neg) and correct.
    for k in range(len(As[0, :])):
        limA = coeff0(q, k)  # limiting value for each entry of eigenvector.
        if _np.sign(limA) != _np.sign(As[0, k]):
            As[0, k] = -As[0, k]
    for k in range(1, len(As[:, 0])):
        for m in range(len(As[0, :])):
            if _np.sign(As[m, n]) != _np.sign(As[m - 1, n]):
                As[m, n] = - As[m, n]
    return As


def coeff0(q, r):
    ''' Limiting value of Fouerier coefficients associated with zeroth
    eigenvalue, for |q| << 1. These are used to check the correct sign of
    the numerically calculate eigenvectors components as a function of q.
    '''
    if r == 0:  # first coefficient
        coeff = 1 / _np.sqrt(2)
    else:
        coeff = ((-q) ** r) / ((4 ** r) * (factorial(r) ** 2))
    return coeff
