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
        As=None,
    ):
        """Cosine elliptic function ce_{2n}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `a`, and the domain.
        """
        # initialize the coefficients (q=0)
        if As is None:
            As = A_coefficients(q, N, type, period)
        vals = {}
        for n in range(N):
            terms = [_np.cos((2 * k) * x) * (As['A' + str(2 * n)][0, k])
                     for k in range(N)]
            vals.update({'ce' + str(2 * n): _np.sum(terms, axis=0)})
            vals.update({'ce' + str(2 * n):
                         vals['ce' + str(2 * n)][_np.newaxis, :]})
            vals.update({'a' + str(2 * n): As['a' + str(2 * n)]})
        for i in range(1, len(q)):
            for n in range(N):
                terms = [_np.cos((2*k)*x) * (As['A' + str(2 * n)][i, k])
                         for k in range(N)]
                ce = _np.sum(terms, axis=0)
                ce = ce[_np.newaxis, :]
                ce = _np.append(vals['ce' + str(2 * n)], ce, axis=0)
                vals.update({'ce' + str(2 * n): ce})
        return vals

    @classmethod
    def ce_odd(
        cls,
        q,
        x,
        N,
        type='odd',
        period='two',
    ):
        """Cosine elliptic function ce_{2n+1}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `a`, and the domain.
        """
        pass

    @classmethod
    def se_even(
        cls,
        q,
        x,
        N,
        type='even',
        period='one',
    ):
        """Sine elliptic function se_{2n+2}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `b`, and the domain.
        """
        pass

    @classmethod
    def se_odd(
        cls,
        q,
        x,
        N,
        type='even',
        period='two',
    ):
        """Sine elliptic function se_{2n+1}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `b`, and the domain.
        """
        pass


def A_coefficients(q, N, type, period):
    vals = {}
    for n in range(N):
        a, A = eig_pairs(matrix_system(q[0], N, type, period))
        a = [a[n]]  # makes a list of the nth eigenvalue
        As = Anorm(A[:, n])
        As = As[_np.newaxis, :]
        for k in range(1, len(q)):
            an, A = eig_pairs(matrix_system(q[k], N))
            a.append(an[n])
            nA = Anorm(A[:, n])
            nAs = nA[_np.newaxis, :]
            As = _np.append(As, nAs, axis=0)
        As = Fcoeffs(As, n)
        vals.update({'a' + str(2 * n): _np.array(a)})
        vals.update({'A' + str(2 * n): As})
    return vals


def solve_system(an, b, q, N, type, period):
    """Inverts the system
        (A - a_n*I)x_{n} = b
    and returns the Fourier coefficients `x_{n}` resulting from the inversion.

    Input:
        q: float(a)
        an: 1d-array. Nth-eigenvalue of Mathieu's equation. a_{n} = a_{n}(q).
            len(a_{n}) = len(q).
        q
        b: 1d-array
    """
    vals = {}
    Ain = matrix_inhom(an[0], q[0], N, type, period)  # matrix to invert
    invA = _np.inv(Ain)
    As = _np.dot(invA * b)  # matrix multiplication on rhs of system
    As = Anorm(As)
    As = As[_np.newaxis, :]
    for k in range(1, len(q)):
        Ain = matrix_inhom(an[k], q[k], N, type, period)
        invA = _np.inv(Ain)
        nAs = _np.dot(invA * b)
        nAs = Anorm(nAs)
        nAs = nAs[_np.newaxis, :]
        As = _np.append(As, nAs, axis=0)
    for n in range(N):
        vals.update({'Ai' + str(2 * n): As[:, n]})
    return vals


def matrix_inhom(an, q, N, type, period):
    """Returns the matrix used to invert when solving an inhomogeneous system,
    associated with a particular eigenvalue. This is
        Ai = A - a_nI
    Input:
        an: value of eigenvalue for fixed q.
        q: Parameter value, either real or (purely) imaginary.
        N: integer, order (size) or matrix. This value is dependent on
            the value of q. For larger values of q, a greater N is needed.
        type: str, either `even` or `odd`.
        period: str, either `one` or `two`. Together with `type`, define the
            eigensolution.
    """
    A = matrix_system(q, N, type, period)
    Ain = A - an * _np.eye(N)
    return Ain

def Fcoeffs(As, n=0, q=0.00001):
    """ Returns the Fourier coefficient of the Mathieu functions for given
    parameter q. Makes sure the coefficients are continuous (in q). Numerical
    routines for estimating eigenvectors might converge in different signs
    for the eigenvectors for different (neighboring) values of q.
        Input:
            As: 2d array. Eigenvector shape(As)=Nq, N, as a function of q and
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
        if n == 0:
            limA = coeff0(q, k)  # limit value for each entry of eig-vector.
        else:
            limA = coeffs(q, k, n)
        if _np.sign(As[0, k]) == 0:
            As[0, k] = limA  # limiting value of coeff for small q
        else:
            if _np.sign(limA) != _np.sign(As[0, k]):
                As[0, k] = -As[0, k]
    for k in range(1, len(As[:, 0])):  # for all values in q
        for m in range(len(As[0, :])):  # iterate through F coeffs
            if _np.sign(As[k, m]) != _np.sign(As[k - 1, m]):
                As[k, m] = -As[k, m]
    for m in range(len(As[0, :])):
        nAs = reflect_coeffs(As[:, m])
        As[:, m] = reflect_coeffs(nAs)  # second crossing of F coeffs
    return As


def reflect_coeffs(A):
    """Fixes the reflection caused on F coefficients when approaching zero,
    in their q-dependence. This typically happens once.
    Input:
        A: 1d array. Fourier coefficient as a function of q.
    Output:
        A: 1d array. Same dimensions and magnitude as original array.
    """
    diffA = _np.gradient(A)  # center order differencing.
    for k in range(1, len(A) - 1):
        if _np.sign(diffA[k + 1]) != _np.sign(diffA[k - 1]):
            if abs(A[k]) < 0.1 and abs(A[k - 1]) > abs(A[k]):
                A[k:] = -A[k:]  # reflect rest of values
                break
    return A


def coeff0(q, r):
    ''' Limiting value of Fourier coefficients associated with zeroth
    eigenvalue, for |q| << 1. These are used to check the correct sign of
    the numerically calculate eigenvectors components as a function of q.
    '''
    if r == 0:  # first coefficient
        coeff = 1 / _np.sqrt(2)
    else:
        coeff = ((-q) ** r) / ((4 ** r) * (factorial(r) ** 2))
    return coeff


def coeffs(q, r, n):
    """Limiting value of Fourier coefficients other than that of zeroth
    eigenvalue. Only for q near zero.
    """
    if r == 0:
        coeff = (1 / (n * factorial((2 * n) - 1))) * (q / 4) ** n
    else:
        if n < r:  # issue here, b/c not always start positive
            nume = ((-1)**(r - n)) * factorial(2 * n)
            denom = factorial(r - n) * factorial(r + n)
            coeff = (nume / denom) * (q / 4) * (r - n)
        elif n == r:
            nume = - (4 * (n ** 2) + 1) / ((4 * (n ** 2) - 1)**2)
            coeff = 1 + (nume * ((q / 4)**2))
        elif n > r:
            nume = factorial(n + r - 1)
            denom = factorial(n - r) * factorial((2 * n) - r)
            coeff = (nume / denom) * (q / 4)**(n - r)
    return coeff


def Anorm(A, type='ce2n'):
    """ Normalization of eigenvectors in accordance to Mathieu functions.
    Default is for that associated with ce_{2n}(q, z).
    Input:
        A: 1d-array. len(A) = N, N being the number of Fourier coefficients.
        type: str, default `ce2n`. Normalization for other functions is
        different.
    Output:
        A: 1d-array. Normalized eigenvector.
    """
    if type is not 'ce2n':
        Astar = _np.conjugate(A)
        norm = _np.sum(A * Astar)
    else:
        A0 = A[0]
        A0star = _np.conjugate(A0)
        A2nstar = _np.conjugate(A[1:])
        norm = 2 * (A0 * A0star) + _np.sum(A[1:] * A2nstar)
    A = A / norm
    return A
