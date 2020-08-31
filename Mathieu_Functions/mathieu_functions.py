"""
defines a class where the Mathieu functions ce and se are defined.
"""
from eig_system import matrix_system, eig_pairs, pseudo_inverse
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
        Ncut=0
    ):
        """Cosine elliptic function ce_{2n}, as a function of parameter
        q (which can be real or purely imaginary), the characteristic number
        `a`, and the domain.
        """
        # initialize the coefficients (q=0)
        if As is None:
            As = A_coefficients(q, N, type, period)
        vals = {}
        if Ncut != 0:
            N = Ncut
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
    if q[0].imag != 0:
        imag = True
    else:
        imag = False
    for n in range(N):
        a, A = eig_pairs(matrix_system(q[0], N, type, period), type, period)
        a = [a[n]]  # makes a list of the nth eigenvalue
        As = Anorm(A[:, n], type, period)
        As = As[_np.newaxis, :]
        for k in range(1, len(q)):
            an, A = eig_pairs(matrix_system(q[k], N), type, period)
            a.append(an[n])
            nA = Anorm(A[:, n], type, period)
            nAs = nA[_np.newaxis, :]
            As = _np.append(As, nAs, axis=0)
        As = Fcoeffs(As, n, flag=imag)
        vals.update({'a' + str(2 * n): _np.array(a)})
        vals.update({'A' + str(2 * n): As})
    return vals


def correction_term(a0, b, eps, q, N, type, period):
    """Calculates the first order correction term
    Input:
        a: 1d-array. It is the eigenvalue of a Mathieu function
        b: 1d-array, inhomogeneous unit-length vector of size 1xN.
        eps: float, real and positive << 1. Perturbation parameter.
        q: 1d-array of float numbers. Real or purely imaginary.
        N: int, the order (size) of the square matrix, whose off-diagonal
        entries contain the parameter q.
        type: str, either `even` or `odd`, referencing the type of Mathieu
            function.
        period: str, either `one` or `two`, refering to pi-periodic or
            two-pi periodic.
    Outout:
        A+: Moore-Penrose matrix that is the pseudo inverse of A-aI.
    """
    A = matrix_system(q[0], N, type, period)
    Aplus = pseudo_inverse(a0[0], q[0], N, type, period)
    X1 = _np.dot(Aplus, b)  # should be orthogonal to x0
    AX1 = _np.dot(A, X1)
    AX1 = AX1 - b
    a1 = [_np.dot(_np.conjugate(X1).T, AX1)[0]]
    if len(q) > 1:  # allows to test for single vals of q.
        X1[0] = X1[0] / _np.sqrt(2)
        X1 = X1[_np.newaxis, :, 0]
        # a1 = [a1]
        for k in range(1, len(q)):
            A = matrix_system(q[k], N, type, period)
            Aplus = pseudo_inverse(a0[k], q[k], N, type, period)
            x1 = _np.dot(Aplus, b)
            Ax1 = _np.dot(A, x1)
            Ax1 = Ax1 - b
            a1.append(_np.dot(_np.conjugate(x1).T, Ax1)[0])
            x1[0] = x1[0] / _np.sqrt(2)
            x1 = x1[_np.newaxis, :, 0]
            X1 = _np.append(X1, x1, axis=0)
    return {'a1': a1, 'X1': X1}


def Fcoeffs(As, n=0, q=0.00001, flag=False):
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
    if flag is True:
        q = q * (1j)
    else:
        for k in range(len(As[0, :])):
            if n == 0:
                limA = coeff0(q, k)  # q~0 value for each entry of eig-vector.
            else:
                limA = coeffs(q, k, n)
            if _np.sign(As[0, k]) == 0:
                As[0, k] = limA  # q ~ 0 value of coeff for small q, n>0
            else:
                if _np.sign(limA) != _np.sign(As[0, k]):
                    As[0, k] = -As[0, k]
        for k in range(1, len(As[:, 0])):  # for all values in q
            for m in range(len(As[0, :])):  # iterate through F coeffs
                if _np.sign(As[k, m]) != _np.sign(As[k - 1, m]):
                    As[k, m] = -As[k, m]
    if flag is False:
        for m in range(len(As[0, :])):
            nAs = reflect_coeffs(As[:, m])
            mAs = reflect_coeffs(nAs)  # second crossing of F coeffs
            As[:, m] = reflect_coeffs(mAs)  # third reflection
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
            power = (q / 4) * (r - n)
            coeff = (nume / denom) * power
            if r >= 8 + n:
                coeff = -coeff
        elif n == r:
            nume = - (4 * (n ** 2) + 1) / ((4 * (n ** 2) - 1)**2)
            coeff = 1 + (nume * ((q / 4)**2))
        elif n > r:
            nume = factorial(n + r - 1)
            denom = factorial(n - r) * factorial((2 * n) - r)
            if n >= 7 and n == r + 1:
                q = 1e-11  # small enough now, when power is order 1
            power = (q / 4)**(n - r)
            coeff = (nume / denom) * power
    return coeff


def Anorm(A, type='even', period='one'):
    """ Normalization of eigenvectors in accordance to Mathieu functions.
    Default is for that associated with ce_{2n}(q, z).
    Input:
        A: 1d-array. len(A) = N, N being the number of Fourier coefficients.
        type: str, default `ce2n`. Normalization for other functions is
        different.
    Output:
        A: 1d-array. Normalized eigenvector.
    """
    if [type, period] == ['even', 'one']:
        A0star = A[0]  # _np.conjugate(A[0])
        Astar = A[1:]  # _np.conjugate(A[1:])
        norm = (2 * (A[0] * A0star)) + _np.sum(A[1:] * Astar)
    else:
        norm = _np.sum(A * _np.conjugate(A))
    A = A / norm
    return A

