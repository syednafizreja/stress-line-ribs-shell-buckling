"""Form-found shell surface z(x, y): force density (ShellBench) plus a bicubic interpolant.

The vault is a height field over the square plan, so everything downstream (meshing, rib
curves, tensor tracing) works in plan and lifts points with ``Surface.z``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RectBivariateSpline
from shellbench.formfind import Vault, form_find

from .constants import G, RHO, RISE, SPAN, T_REF


@dataclass(frozen=True)
class Surface:
    """Bicubic surface through the force-density nodes."""

    vault: Vault
    spline: RectBivariateSpline  # evaluated as spline(y, x): rows are y
    span: float

    def z(self, x, y):
        return self.spline.ev(np.asarray(y, float), np.asarray(x, float))

    def grad(self, x, y):
        """Return (dz/dx, dz/dy)."""
        y, x = np.asarray(y, float), np.asarray(x, float)
        return self.spline.ev(y, x, dy=1), self.spline.ev(y, x, dx=1)

    def normal(self, x, y):
        """Unit upward normal, shape (..., 3)."""
        zx, zy = self.grad(x, y)
        n = np.stack([-zx, -zy, np.ones_like(zx)], axis=-1)
        return n / np.linalg.norm(n, axis=-1, keepdims=True)


def shell_surface(span: float = SPAN, rise: float = RISE, n: int = 60) -> Surface:
    """Form-find under the self-weight of the reference shell and return the surface.

    The shape does not depend on the load magnitude: the force density is rescaled so that the
    crown reaches ``rise``.
    """
    v = form_find(span, span, rise, w=RHO * G * T_REF, nx=n, ny=n)
    spline = RectBivariateSpline(v.y[:, 0], v.x[0, :], v.z, kx=3, ky=3)
    return Surface(v, spline, span)


def _corner_fdm(span: float, n: int, edge_ratio: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fixed-plan force density, unit q, only the 4 corners fixed; returns x, y, z (hanging form, z >= 0).

    Boundary branches carry ``edge_ratio`` x the interior force density (the edge arches).
    Tributary nodal load is w*dx*dy with w = 1 (z is rescaled afterwards).
    """
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import spsolve

    m = n + 1
    xs = np.linspace(0.0, span, m)
    x, y = np.meshgrid(xs, xs)
    idx = lambda i, j: i * m + j  # noqa: E731  (row i = y, col j = x)
    fixed = {idx(0, 0), idx(0, n), idx(n, 0), idx(n, n)}
    d = span / n
    a = lil_matrix((m * m, m * m))
    p = np.zeros(m * m)
    for i in range(m):
        for j in range(m):
            k = idx(i, j)
            trib = d * d * (0.5 if i in (0, n) else 1.0) * (0.5 if j in (0, n) else 1.0)
            p[k] = trib
            for di, dj in ((0, 1), (1, 0)):
                ii, jj = i + di, j + dj
                if ii > n or jj > n:
                    continue
                on_edge = (di == 0 and i in (0, n)) or (dj == 0 and j in (0, n))
                q = (edge_ratio if on_edge else 1.0) / d
                k2 = idx(ii, jj)
                a[k, k] += q
                a[k2, k2] += q
                a[k, k2] -= q
                a[k2, k] -= q
    free = np.array([k for k in range(m * m) if k not in fixed])
    z = np.zeros(m * m)
    z[free] = spsolve(a.tocsr()[free][:, free], p[free])  # sum q (z_i - z_j) = p_i (inverted hanging form)
    return x, y, z.reshape(m, m)


def corner_surface(span: float = SPAN, rise: float = RISE, n: int = 60, edge_rise_frac: float = 0.5) -> Surface:
    """Corner-supported funicular shell (scope section 10).

    Fixed-plan force density with only the corners fixed. The boundary force density is chosen so
    the edge-arch rise is ``edge_rise_frac`` x the crown rise; the free edges then need an edge member
    to carry their in-plane thrust (modelled as an edge beam in ``fe``).
    """
    from scipy.optimize import brentq

    def frac(r: float) -> float:
        _, _, z = _corner_fdm(span, n, r)
        return z[0, n // 2] / z[n // 2, n // 2] - edge_rise_frac

    r = brentq(frac, 0.5, 500.0, xtol=1e-4)
    x, y, z = _corner_fdm(span, n, r)
    z *= rise / z.max()
    d = span / n
    v = Vault(x, y, z, q_x=float("nan"), q_y=float("nan"), dx=d, dy=d, w=RHO * G * T_REF)
    spline = RectBivariateSpline(y[:, 0], x[0, :], z, kx=3, ky=3)
    return Surface(v, spline, span)
