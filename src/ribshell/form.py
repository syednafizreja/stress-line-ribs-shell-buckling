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
