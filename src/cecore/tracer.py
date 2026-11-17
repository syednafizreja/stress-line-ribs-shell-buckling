"""Principal-direction (tensor) line tracer for 2-D symmetric tensor fields on a rectangle.

Used by P2 (rib layouts from membrane-force and moment fields) and P1 Part B.

Method
------
* The field is given at scattered plan points (element centroids) and resampled on a regular grid
  (linear interpolation, nearest-neighbour fill near the boundary), then read with bilinear
  interpolation of the three independent components.
* Lines follow the eigenvector of the larger (``which="max"``) or smaller (``"min"``) eigenvalue.
  Tensor lines have no orientation, so at every RK4 stage the eigenvector sign is flipped to agree
  with the previous direction.
* A line stops at the boundary (its last point is moved onto it), near degenerate points where the
  anisotropy |l1 - l2| / (|l1| + |l2|) falls below ``degenerate_tol``, when it comes closer than
  ``d_test`` to another line of the same family, where the tensor magnitude (Frobenius norm) drops
  below ``magnitude_tol`` x its field maximum (directions there are noise and carry no load worth a
  rib), when it closes on itself, or after ``max_len``.
* Seeding follows Jobard & Lefer (1997): the first seed is the most anisotropic grid point; new seed
  candidates sit at +-d_sep normal to every accepted line, and a candidate is used only when it is at
  least d_sep from every line of the family.
* Output polylines are resampled at ``ds_out``, lightly smoothed with fixed end points, and end
  points closer than ``snap`` to the boundary are moved onto it.

Reference: Jobard, B. & Lefer, W. (1997) Creating evenly-spaced streamlines of arbitrary density.
Visualization in Scientific Computing '97, 43-55. doi:10.1007/978-3-7091-6876-9_5
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, RegularGridInterpolator
from scipy.spatial import cKDTree


@dataclass
class TensorField:
    """Bilinear field of 2 x 2 symmetric tensors on the plan rectangle [x0, x1] x [y0, y1]."""

    x0: float
    x1: float
    y0: float
    y1: float
    _f: RegularGridInterpolator

    @classmethod
    def from_points(cls, xy: np.ndarray, tensors: np.ndarray, bounds, n: int = 121) -> "TensorField":
        x0, x1, y0, y1 = bounds
        comp = np.stack([tensors[:, 0, 0], tensors[:, 0, 1], tensors[:, 1, 1]], axis=-1)
        gx, gy = np.linspace(x0, x1, n), np.linspace(y0, y1, n)
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        q = np.column_stack([X.ravel(), Y.ravel()])
        vals = LinearNDInterpolator(xy, comp)(q)
        bad = np.isnan(vals[:, 0])
        if bad.any():
            vals[bad] = NearestNDInterpolator(xy, comp)(q[bad])
        return cls(x0, x1, y0, y1, RegularGridInterpolator((gx, gy), vals.reshape(n, n, 3)))

    def symmetrized(self, group: str | None) -> "TensorField":
        """Group-average of the field: T_s(p) = mean_g R_g^T T(c + R_g (p - c)) R_g (square domain)."""
        if not group:
            return self
        gx, gy = self._f.grid
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        c = np.array([0.5 * (self.x0 + self.x1), 0.5 * (self.y0 + self.y1)])
        p = np.column_stack([X.ravel(), Y.ravel()]) - c
        acc = np.zeros((p.shape[0], 2, 2))
        for R in symmetry_ops(group):
            a, b, cc = self._f(np.clip(p @ R.T + c, [self.x0, self.y0], [self.x1, self.y1])).T
            T = np.stack([np.stack([a, b], -1), np.stack([b, cc], -1)], -2)
            acc += np.einsum("ji,njk,kl->nil", R, T, R)
        acc /= len(symmetry_ops(group))
        vals = np.stack([acc[:, 0, 0], acc[:, 0, 1], acc[:, 1, 1]], -1).reshape(len(gx), len(gy), 3)
        return TensorField(self.x0, self.x1, self.y0, self.y1, RegularGridInterpolator((gx, gy), vals))

    def inside(self, p) -> bool:
        return self.x0 <= p[0] <= self.x1 and self.y0 <= p[1] <= self.y1

    def eig(self, p):
        """Eigenvalues (ascending) and eigenvectors (columns) at plan point p."""
        a, b, c = self._f(np.clip([[p[0], p[1]]], [self.x0, self.y0], [self.x1, self.y1]))[0]
        return np.linalg.eigh(np.array([[a, b], [b, c]]))

    def magnitude(self, p) -> float:
        w, _ = self.eig(p)
        return float(np.hypot(w[0], w[1]))

    def max_magnitude(self, n: int = 61) -> float:
        g = [(x, y) for x in np.linspace(self.x0, self.x1, n) for y in np.linspace(self.y0, self.y1, n)]
        return max(self.magnitude(p) for p in g)

    def anisotropy(self, p) -> float:
        w, _ = self.eig(p)
        return abs(w[1] - w[0]) / (abs(w[0]) + abs(w[1]) + 1e-30)

    def direction(self, p, which: str, prev=None):
        _, v = self.eig(p)
        d = v[:, 1] if which == "max" else v[:, 0]
        if prev is not None and np.dot(d, prev) < 0:
            d = -d
        return d


def symmetry_ops(group: str | None) -> list[np.ndarray]:
    """2 x 2 orthogonal maps about the domain centre: "D4" (square), "mirror_y" (v -> -v), "mirror_x"."""
    I = np.eye(2)
    if not group:
        return [I]
    if group == "mirror_y":
        return [I, np.diag([1.0, -1.0])]
    if group == "mirror_x":
        return [I, np.diag([-1.0, 1.0])]
    if group == "D4":
        S = np.array([[0.0, 1.0], [1.0, 0.0]])
        flips = [np.diag(d) for d in ((1, 1), (-1, 1), (1, -1), (-1, -1))]
        return flips + [F @ S for F in flips]
    raise ValueError(f"unknown symmetry group {group!r}")


def _orbit(line: np.ndarray, f: TensorField, group: str | None, tol: float) -> list[np.ndarray]:
    """Distinct images of ``line`` under the group (an image that retraces a kept one is dropped)."""
    c = np.array([0.5 * (f.x0 + f.x1), 0.5 * (f.y0 + f.y1)])
    out: list[np.ndarray] = []
    for R in symmetry_ops(group):
        img = (line - c) @ R.T + c
        if all(np.median(cKDTree(o).query(img)[0]) > tol for o in out):
            out.append(img)
    return out


def _trim_to_orbit(line: np.ndarray, f: TensorField, group: str | None, tol: float, d_test: float
                   ) -> np.ndarray:
    """Longest run of ``line`` that stays d_test away from its own non-coincident images.

    The orbit of the trimmed line is then mutually separated: dist(g1 L', g2 L') = dist(L', g L')
    for g = g1^-1 g2, and L' keeps d_test from every g L that does not retrace L.
    """
    c = np.array([0.5 * (f.x0 + f.x1), 0.5 * (f.y0 + f.y1)])
    own = cKDTree(line)
    imgs = []
    for R in symmetry_ops(group):
        img = (line - c) @ R.T + c
        if np.median(own.query(img)[0]) > tol:      # not a symmetry of this line
            imgs.append(img)
    if not imgs:
        return line
    ok = cKDTree(np.vstack(imgs)).query(line)[0] >= d_test
    best, start, run = (0, 0), None, []
    for i, flag in enumerate(np.append(ok, False)):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return line[best[0]:best[1]]


def _boundary_hit(p, q, f: TensorField):
    """Point where segment p -> q leaves the rectangle."""
    t = 1.0
    for k, lo, hi in ((0, f.x0, f.x1), (1, f.y0, f.y1)):
        dq = q[k] - p[k]
        if q[k] < lo:
            t = min(t, (lo - p[k]) / dq)
        elif q[k] > hi:
            t = min(t, (hi - p[k]) / dq)
    return p + max(0.0, t) * (q - p)


def _integrate(f, seed, d0, which, step, d_test, trees, degenerate_tol, max_len, mag_min=0.0):
    pts, d, length = [np.asarray(seed, float)], np.asarray(d0, float), 0.0
    while length < max_len:
        p = pts[-1]
        k1 = f.direction(p, which, d)
        k2 = f.direction(p + 0.5 * step * k1, which, k1)
        k3 = f.direction(p + 0.5 * step * k2, which, k2)
        k4 = f.direction(p + step * k3, which, k3)
        dn = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        dn /= np.linalg.norm(dn)
        q = p + step * dn
        if not f.inside(q):
            pts.append(_boundary_hit(p, q, f))
            break
        if f.anisotropy(q) < degenerate_tol or (mag_min and f.magnitude(q) < mag_min):
            break
        if any(t is not None and t.query(q)[0] < d_test for t in trees):
            break
        if len(pts) > 20 and np.linalg.norm(q - pts[0]) < 0.5 * step:  # closed loop
            break
        pts.append(q)
        d, length = dn, length + step
    return pts


def _trace_both(f, seed, which, step, d_test, trees, tol, max_len, mag_min=0.0):
    d0 = f.direction(seed, which)
    fwd = _integrate(f, seed, d0, which, step, d_test, trees, tol, max_len, mag_min)
    bwd = _integrate(f, seed, -d0, which, step, d_test, trees, tol, max_len, mag_min)
    return np.array(bwd[::-1] + fwd[1:])


def _resample(line: np.ndarray, ds: float) -> np.ndarray:
    seg = np.linalg.norm(np.diff(line, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    n = max(2, int(round(s[-1] / ds)) + 1)
    t = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(t, s, line[:, 0]), np.interp(t, s, line[:, 1])])


def _smooth(line: np.ndarray, iters: int = 3) -> np.ndarray:
    out = line.copy()
    for _ in range(iters):
        out[1:-1] = 0.25 * out[:-2] + 0.5 * out[1:-1] + 0.25 * out[2:]
    return out


def _snap(p, f: TensorField, snap: float):
    p = p.copy()
    for k, lo, hi in ((0, f.x0, f.x1), (1, f.y0, f.y1)):
        if p[k] - lo < snap:
            p[k] = lo
        elif hi - p[k] < snap:
            p[k] = hi
    return p


def trace_lines(points_xy: np.ndarray, tensors: np.ndarray, bounds, d_sep: float = 1.5,
                step: float = 0.05, degenerate_tol: float = 0.05, which: str = "max",
                d_test_frac: float = 0.5, min_len_frac: float = 1.0, ds_out: float = 0.5,
                snap: float = 0.3, max_len: float = 100.0, field: TensorField | None = None,
                magnitude_tol: float = 0.0, symmetry: str | None = None
                ) -> list[np.ndarray]:
    """Evenly spaced tensor lines of one family.

    Args:
        points_xy: (n, 2) plan points where ``tensors`` are given.
        tensors: (n, 2, 2) symmetric tensors in plan components.
        bounds: (x0, x1, y0, y1).
        d_sep: line separation (Jobard-Lefer); d_test = d_test_frac * d_sep.
        which: "max" or "min" eigenvalue family.
        min_len_frac: lines shorter than min_len_frac * d_sep are dropped.
        magnitude_tol: lines stop where |T| < magnitude_tol * max |T| (0 disables).
        symmetry: None, "D4", "mirror_y" or "mirror_x". The field is group-averaged (removes
            mesh-level asymmetry) and every accepted line is added with its whole orbit, so the
            layout has the symmetry of the structure and load. Each line is first trimmed to its
            longest run that stays d_test from its own images, so the orbit stays separated.

    Returns:
        Plan polylines (k x 2 arrays), resampled at about ``ds_out``.
    """
    f = (field or TensorField.from_points(points_xy, tensors, bounds)).symmetrized(symmetry)
    d_test = d_test_frac * d_sep
    mag_min = magnitude_tol * f.max_magnitude() if magnitude_tol else 0.0
    # first seed: most anisotropic point of a coarse grid
    gx = np.linspace(f.x0, f.x1, 25)[1:-1]
    cand = [np.array([x, y]) for x in gx for y in gx]
    seed0 = max(cand, key=lambda p: f.anisotropy(p) * f.magnitude(p))
    lines: list[np.ndarray] = []
    trees: list[cKDTree | None] = []
    queue = [seed0]
    while queue:
        s = queue.pop(0)
        if not f.inside(s) or f.anisotropy(s) < degenerate_tol or (mag_min and f.magnitude(s) < mag_min):
            continue
        if trees and min(t.query(s)[0] for t in trees) < d_sep * 0.999:
            continue
        line = _trace_both(f, s, which, step, d_test, trees, degenerate_tol, max_len, mag_min)
        if len(line) < 2 or np.sum(np.linalg.norm(np.diff(line, axis=0), axis=1)) < min_len_frac * d_sep:
            continue
        if symmetry:
            line = _trim_to_orbit(line, f, symmetry, 0.1 * d_sep, d_test)
            if len(line) < 2 or np.sum(np.linalg.norm(np.diff(line, axis=0), axis=1)) < min_len_frac * d_sep:
                continue
        orbit = _orbit(line, f, symmetry, 0.1 * d_sep)
        for ln in orbit:
            lines.append(ln)
            trees.append(cKDTree(ln))
            # new seed candidates at +-d_sep normal to the line, every d_sep / 2 along it
            for p in _resample(ln, 0.5 * d_sep):
                nrm = f.direction(p, "min" if which == "max" else "max")
                queue += [p + d_sep * nrm, p - d_sep * nrm]
    out = []
    for ln in lines:
        r = _smooth(_resample(ln, ds_out))
        r[0], r[-1] = _snap(r[0], f, snap), _snap(r[-1], f, snap)
        out.append(r)
    return out


def line_stats(lines: list[np.ndarray]) -> dict[str, float]:
    """Total length and minimum radius of curvature (three-point circle) over all lines."""
    total, rmin = 0.0, np.inf
    for ln in lines:
        total += float(np.sum(np.linalg.norm(np.diff(ln, axis=0), axis=1)))
        for a, b, c in zip(ln[:-2], ln[1:-1], ln[2:]):
            ab, bc, ca = np.linalg.norm(b - a), np.linalg.norm(c - b), np.linalg.norm(a - c)
            area2 = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            if area2 > 1e-12:
                rmin = min(rmin, ab * bc * ca / (2.0 * area2))
    return {"length": total, "min_radius": rmin}
