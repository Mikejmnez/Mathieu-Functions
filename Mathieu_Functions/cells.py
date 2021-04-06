"""
defines a class where the cellular flow solution is calculated
"""
import numpy as _np
import xarray as _xr

class cellular_flow:

    def __init__(self, q, N):
        self._q = q
        self._N = N

    @classmethod
    def cells1(
        cls,
        Q,
        k,
        x,
        y,
        N,
        As,
        lc=0,
    ):
        """Calculates the O(epsilon) solution to 2d advection diffusion
        equation.
        returns a dictionary containint the all N O(e) functions that
        constitute the next order solution (in their n-sum)
        """

        #  iterates for several values of k
        ds_r, ds_i = cell1_iterator(As, k, y, x, Q, N, lc)
        return ds_r, ds_i

    @classmethod
    def cells2(
        cls,
        Q,
        k,
        x,
        y,
        N,
        As,
        lc=0,
    ):
        """Calculates the O(epsilon) solution to 2d advection diffusion
        equation.
        returns a dictionary containint the all N O(e) functions that
        constitute the next order solution (in their n-sum)
        """

        #  iterates for several values of k
        ds_r, ds_i = cell2_iterator(As, k, y, x, Q, N, lc)
        return ds_r, ds_i


def cell1_iterator(As, k, y, x, Q, N, lc=0):
    if type(Q) == _np.complex128:
        Q = [Q]
        nq = -1
    else:
        nq = 0
    vars_r = []
    vars_i = []
    coords = {"q": Q,
              "y": 2 * y[:, 0],
              "x": x[0, :]}
    for n in range(N):
        SUM = []
        for r in range(N):
            A2r = As["A" + str(2 * n)][nq, r]
            f1 = (k + r) / 4
            f2 = (k - r) / 4
            r11 = f1 * _np.cos(2 * (lc - r) * y)
            r12 = f2 * _np.cos(2 * (lc + r) * y)
            r1 = A2r * (r11 + r12) * _np.exp(1j * (k + lc) * x)
            r21 = f2 * _np.cos(2 * (lc - r) * y)
            r22 = f1 * _np.cos(2 * (lc + r) * y)
            r2 = A2r * (r21 + r22) * _np.exp(1j * (k - lc) * x)
            SUM.append(r1 + r2)
        nterm = sum(SUM)
        data_r = _xr.DataArray(coords=coords, dims=["q", "y", "x"])
        data_i = _xr.DataArray(coords=coords, dims=["q", "y", "x"])
        data_r.sel(q=Q[0])[:] = nterm.real
        data_i.sel(q=Q[0])[:] = nterm.imag
        vars_r.append(_xr.Dataset({'theta1_' + str(2 * n): data_r}))
        vars_i.append(_xr.Dataset({'theta1_' + str(2 * n): data_i}))
    ds_r = _xr.merge(vars_r)  # dataset with real part of the sum
    ds_i = _xr.merge(vars_i)  # dataset with imag part of the sum
    for i in range(1, len(Q)):  # only loops if Q has more that one element
        for n in range(N):
            SUM = []
            for r in range(N):
                A2r = As['A' + str(2 * n)][i, r]  # n=0
                f1 = (k + r) / 4
                f2 = (k - r) / 4
                r11 = f1 * _np.cos(2 * (lc - r) * y)
                r12 = f2 * _np.cos(2 * (lc + r) * y)
                r1 = A2r * (r11 + r12) * _np.exp(1j * (k + lc) * x)
                r21 = f2 * _np.cos(2 * (lc - r) * y)
                r22 = f1 * _np.cos(2 * (lc + r) * y)
                r2 = A2r * (r21 + r22) * _np.exp(1j * (k - lc) * x)
                SUM.append(r1 + r2)
            nterm = sum(SUM)
            ds_r['theta1_' + str(2 * n)].sel(q=Q[i].imag)[:] = nterm.real
            ds_i['theta1_' + str(2 * n)].sel(q=Q[i].imag)[:] = nterm.imag
    return ds_r, ds_i


def cell2_iterator(As, k, y, x, Q, N, lc=0):
    if type(Q) == _np.complex128:
        Q = [Q]
        nq = -1
    else:
        nq = 0
    vars_r = []
    vars_i = []
    coords = {"q": Q,
              "y": 2 * y[:, 0],
              "x": x[0, :]}
    for n in range(N):
        SUM = []
        for r in range(N):
            A2r = As["A" + str(2 * n)][nq, r]
            f1_m1 = (k**2 - 2 * k * lc - r**2) / 8
            f2_m1 = (k - r)**2 / 16
            f3_m1 = (k + r)**2 / 16
            r1_m1 = f1_m1 * _np.cos(2 * r * y)
            r2_m1 = f2_m1 * _np.cos(2 * (r - 2 * lc) * y)
            r3_m1 = f3_m1 * _np.cos(2 * (r + 2 * lc) * y)
            rm1 = -A2r * (r1_m1 + r2_m1 + r3_m1) * _np.exp(1j * (k - 2*lc)*x)
            f1_0 = (k**2 + r**2) / 4
            f2_0 = (k**2 + 2 * r * lc - r**2) / 8
            f3_0 = (k**2 - 2 * r * lc - r**2) / 8
            r1_0 = f1_0 * _np.cos(2 * r * y)
            r2_0 = f2_0 * _np.cos(2 * (r - 2 * lc) * y)
            r3_0 = f3_0 * _np.cos(2 * (r + 2 * lc) * y)
            r0 = -A2r * (r1_0 + r2_0 + r3_0) * _np.exp(1j * k * x)
            f1_1 = (k**2 + 2 * k * lc - r**2) / 8
            f2_1 = (k + r)**2 / 16
            f3_1 = (k - r)**2 / 16
            r1_1 = f1_1 * _np.cos(2 * r * y)
            r2_1 = f2_1 * _np.cos(2 * (r - 2 * lc) * y)
            r3_1 = f3_1 * _np.cos(2 * (r + 2 * lc) * y)
            r1 = -A2r * (r1_1 + r2_1 + r3_1) * _np.exp(1j * (k + 2 * lc) * x)
            SUM.append(rm1 + r0 + r1)
        nterm = sum(SUM)
        data_r = _xr.DataArray(coords=coords, dims=["q", "y", "x"])
        data_i = _xr.DataArray(coords=coords, dims=["q", "y", "x"])
        data_r.sel(q=Q[0])[:] = nterm.real
        data_i.sel(q=Q[0])[:] = nterm.imag
        vars_r.append(_xr.Dataset({'theta2_' + str(2 * n): data_r}))
        vars_i.append(_xr.Dataset({'theta2_' + str(2 * n): data_i}))
    ds_r = _xr.merge(vars_r)  # dataset with real part of the sum
    ds_i = _xr.merge(vars_i)  # dataset with imag part of the sum
    # for i in range(1, len(Q)):  # only loops if Q has more that one element
    #     for n in range(N):
    #         SUM = []
    #         for r in range(N):
    #             A2r = As['A' + str(2 * n)][i, r]  # n=0
    #             f1 = (k + r) / 4
    #             f2 = (k - r) / 4
    #             r11 = f1 * _np.cos(2 * (lc - r) * y)
    #             r12 = f2 * _np.cos(2 * (lc + r) * y)
    #             r1 = A2r * (r11 + r12) * _np.exp(1j * (k + lc) * x)
    #             r21 = f2 * _np.cos(2 * (lc - r) * y)
    #             r22 = f1 * _np.cos(2 * (lc + r) * y)
    #             r2 = A2r * (r21 + r22) * _np.exp(1j * (k - lc) * x)
    #             SUM.append(r1 + r2)
    #         nterm = sum(SUM)
    #         ds_r['theta1_' + str(2 * n)].sel(q=Q[i].imag)[:] = nterm.real
    #         ds_i['theta1_' + str(2 * n)].sel(q=Q[i].imag)[:] = nterm.imag
    return ds_r, ds_i
