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
            As = Fcoeffs(As)
            vals.update({'a' + str(2 * n): _np.array(a)})
            vals.update({'A' + str(2 * n): As})
            # ce = As[n, 0]
            # for k in range(1, N):
            #     ce = ce + As[n, k] * _np.cos((2 * k) * x)
            # vals.update({'ce' + str(2 * n): ce})
        return vals


def Fcoeffs(As):
    """ Returns the Fourier coefficient of the Mathieu functions for given
    parameter q.
        Output:
            Dictionary.
    """
    for k in range(1, len(As[:, 0])):
        for n in range(len(As[0, :])):
            if _np.sign(As[k, n]) != _np.sign(As[k - 1, n]):
                As[k, n] = - As[k, n]
    return As


def sign_check(A0, A1):
    """
    Makes sure Fourier coefficients are continuous. Numerical rutines for
    estimating eigenvectors might converge in different signs for the
    eigenvectors for different parameter q (real or purely imaginary).
    Input:
        Dictionary .
    """
    for n in range(1, len(A)):
        if _np.sign(A[n].real) != _np.sign(A[n - 1]):
            A[n] = - A[n]
    return A



def limit_coeff(n, A, q):
    """
    Checks that the sign of the eigenvector components matches the sign of the
    Fourier coefficients associated with each eigenvalue with index n, and given (real or purely imaginary) parameter q.
    Input:
        n: int. position (index) of eigenvalue.
        A: 1-d array. Eigenvector associated with eigenvalue a_{n}.
        q: float, purely real or purely imaginary
    """
    pass


def coeff0(q, r):
    ''' Fouerier coefficients associated with zeroth-eigenvalue, valid only
    when |q| << 1. These are used to check the correct sign of the numerically
    calculate eigenvectors components as a function of q.
    '''
    if r == 0:
        coeff = 1 / np.sqrt(2)
    else:
        coeff = np.sqrt(2)*((-q)**r)/((4**r)*(factorial(r)**2))
    return coeff









