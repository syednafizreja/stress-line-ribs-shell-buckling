"""Variant E: zonal shell thickness at fixed volume.

The plan square is split into nz x nz equal zones; every element takes the thickness of the zone
that holds its centroid. Design vectors are mapped to zone thicknesses either one-to-one (free,
nz^2 variables) or through the D4 orbit of each zone (``d4``: 6 variables for nz = 5).

Equal volume: the zone thicknesses are projected onto
    sum_z A_z t_z = V,   T_MIN <= t_z <= T_MAX
by t_z = clip(s * u_z) with the scalar s found by bisection (monotone in s). This keeps the
search unconstrained (CMA-ES works in u), and every evaluated design has exactly the volume V.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .constants import SPAN, T_MAX, T_MIN
from .mesh import ShellMesh


def zone_of_elements(mesh: ShellMesh, nz: int, span: float = SPAN) -> dict[int, int]:
    cen = mesh.centroids()
    out = {}
    for e, c in cen.items():
        i = min(int(c[0] / span * nz), nz - 1)
        j = min(int(c[1] / span * nz), nz - 1)
        out[e] = i * nz + j
    return out


def d4_orbits(nz: int) -> np.ndarray:
    """zone -> orbit index under the square's symmetry group (folded, unordered (i, j))."""
    keys, idx = {}, np.empty(nz * nz, int)
    for i in range(nz):
        for j in range(nz):
            a, b = min(i, nz - 1 - i), min(j, nz - 1 - j)
            k = (min(a, b), max(a, b))
            idx[i * nz + j] = keys.setdefault(k, len(keys))
    return idx


class ZonalDesign:
    def __init__(self, mesh: ShellMesh, nz: int, volume: float, symmetry: str | None = None):
        self.mesh, self.nz, self.volume = mesh, nz, volume
        self.zone = zone_of_elements(mesh, nz)
        areas = mesh.element_areas()
        self.zone_area = np.zeros(nz * nz)
        for e, z in self.zone.items():
            self.zone_area[z] += areas[e]
        self.map = d4_orbits(nz) if symmetry == "d4" else np.arange(nz * nz)
        self.n_vars = int(self.map.max()) + 1
        if not (T_MIN * self.zone_area.sum() <= volume <= T_MAX * self.zone_area.sum()):
            raise ValueError("volume outside the thickness bounds")

    def zone_thickness(self, u: np.ndarray) -> np.ndarray:
        """Design vector u (> 0, any scale) -> zone thicknesses at exactly the target volume."""
        w = np.maximum(np.asarray(u, float)[self.map], 1e-6)
        vol = lambda s: float(self.zone_area @ np.clip(s * w, T_MIN, T_MAX)) - self.volume  # noqa: E731
        s = brentq(vol, 1e-9, T_MAX / w.min() * 1.0001, xtol=1e-14)
        return np.clip(s * w, T_MIN, T_MAX)

    def element_thickness(self, tz: np.ndarray) -> dict[int, float]:
        return {e: float(tz[z]) for e, z in self.zone.items()}

    def lift(self, u_sym: np.ndarray, other: "ZonalDesign") -> np.ndarray:
        """Express a design of ``other`` (e.g. d4) as a design vector of this (free) parametrisation."""
        return np.asarray(u_sym, float)[other.map]
