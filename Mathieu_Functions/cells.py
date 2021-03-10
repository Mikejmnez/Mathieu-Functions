"""
defines a class where the cellular flow solution is calculated
"""

from eig_system import matrix_system, eig_pairs
import numpy as _np
from mathieu_functions import A_coefficients
from mathieu_functions import mathieu_functions as mfs
import xarray as _xr

class cellular_flow:

    def __init__(self, q, N):
        self._q = q
        self._N = N

    @classmethod
    def cells1(
        cls,
        lc,
        k,
        Pe,
        x,
        y,
        N,
        As,
    ):
        """Calculates the O(epsilon) solution to 2d advection diffusion
        equation.
        returns a dictionary containint the all N O(e) functions that
        constitute the next order solution (in their n-sum)
        """

        #  add many cellular flows and, for now, only a single value of q
        Q = 2 * k * Pe
        coords = {"q": Q,
                  "y": 2 * y,
                  "x": x}
        vars_r = []
        vars_i = []
        for n in range(N):
            SUM = []
            for r in range(N):
                A2r = As["A" + str(2 * n)][0, r]
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
        return ds_r, ds_i
