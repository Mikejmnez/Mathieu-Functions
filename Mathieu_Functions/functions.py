from numpy import linalg as _LA
import numpy as _np


def matrix_system(q, N, type='even', period='one'):
    ''' Creates tridianogal matrix system of order NxN, associated with
    each of the four classes of simply-periodic functions.

    Input:
        q: parameter, real or purely imaginary.
        N: Size of the matrix, and thus the order of the highest harmonic in
            the trigonometric series that defines each Mathieu function.
        type: str, `even` or `odd`.
        period: str, `one` or 'two'. If `one`, function is assumed to be
            pi-periodic. If `two` function is taken to be `2pi`-periodic.
    Outout:
        A: ndarray, the square matrix associated with each of the four types
            of simply-periodic Mathieu-functions.
    '''
    if type is 'even' and period is 'one':
        d = [(2. * r) ** 2 for r in range(N)]
        e = q * _np.ones(N - 1)
        A = _np.diag(d) + _np.diag(e, k=-1) + _np.diag(e, k=1)
        A[0, 1] = _np.sqrt(2) * A[0, 1]
        A[1, 0] = A[0, 1]

    return A


def eig_pairs(A, type='real'):
    ''' Calculates the characteristic value (eigenvalue) and the Fourier
    coefficients associated with the Mathieu function. Both the eigenvalues
    and Fourier coefficients are given in ascending order.
    '''
    N = len(A[:, 0])  # size of square matrix.
    w, V = _LA.eig(A)  # calculates the eigenvalue
    V[0, :] = V[0, :] / np.sqrt(2)
    Fcoeff = []
    for n in range(N):
        Fcoeff.append(abs(V[0, n]) / _np.sqrt(2))

    return w, Fcoeff


