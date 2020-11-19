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
    if type is 'even':
        if period is 'one':  # ce_{2n}
            d = [(2. * r) ** 2 for r in range(N)]
            e = q * _np.ones(N - 1)
            A = _np.diag(d) + _np.diag(e, k=-1) + _np.diag(e, k=1)
            A[0, 1] = _np.sqrt(2) * A[0, 1]
            A[1, 0] = A[0, 1]
        elif period is 'two':  # se_{2n+2}
            pass
    elif type is 'odd':
        if period is 'one':  # se_{2n+1}
            pass
        elif period is 'two':  # ce_{2n+1}
            d = [((2. * r) + 1)**2 for r in range(N)]
            e = q * _np.ones(N - 1)
            A = _np.diag(d) + _np.diag(e, k=-1) + _np.diag(e, k=1)
    return A


def eig_pairs(A, type='even', period='one'):
    ''' Calculates the characteristic value (eigenvalue) and the Fourier
    coefficients associated with the Mathieu function. Both the eigenvalues
    and Fourier coefficients are given in ascending order.
    '''
    w, V = _LA.eig(A)  # calculates the eigenvalue
    if [type, period] == ['even', 'one']:
        V[0, :] = V[0, :] / _np.sqrt(2)  # remove factor on all first entries
    #  Sort the eigenvalues and in accordance, re-arrange eigenvectors
    ord_w, V = order_check(w, V)
    return ord_w, V


def order_check(a, v):
    """ Check the ordering of the eigenvalue array, from smaller to larger. If
    true, return a unchanged. Ordering also matters if a is complex. If a is
    complex, ordering again is first set according to real(a). If two
    eigenvalues are complex conjugates, then ordering is in accordance to the
    sign of complex(a). Negative sign is first.
    """
    if a.imag.any() == 0:
        ordered_a = a
        nv = v
    else:
        Ind = _np.argsort(_np.round(a, 1))  # sorting through 5 decimals
        ordered_a = a[Ind]
        nv = 0 * _np.copy(v)
        for k in range(len(Ind)):
            nv[:, k] = v[:, Ind[k]]
    return ordered_a, nv
