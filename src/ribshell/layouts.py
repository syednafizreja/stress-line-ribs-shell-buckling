"""Rib layouts B, C, D and the equal-volume rib depth.

* B: orthogonal grid at about 1.5 m spacing (no mechanics).
* C: both families of the gravity membrane-force lines (principal directions of n under G).
* D: both families of the principal-moment lines of the seismic envelope tensor
  M_env = sum over (+x, -x, +y, -y) of sum_i |m_i| v_i v_i^T  (brief 5.3).
  D_x (tracing from the +x case only) is kept for the comparison the brief asks for.

Rib depth: the base shell is 50 mm and ribs are 150 mm wide. With A_s the shell area, L the rib
length on the surface and k the number of rib crossings (where two ribs share a 150 x 150 mm
footprint), equal volume V gives

    V = A_s t_b + w d (L - k w)   =>   d = (V - A_s t_b) / (w (L - k w)).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from cecore.tracer import line_stats, trace_lines
from shapely.geometry import LineString

from .constants import RIB_WIDTH, SPAN, T_BASE_RIBBED
from .fe import Model, Rib
from .fields import Resultants, abs_tensor, envelope
from .mesh import build_mesh

D_SEP = 1.5
EDGE_CLEAR = 0.5   # [m]


@dataclass
class Layout:
    name: str
    lines: list[np.ndarray]   # plan polylines


def grid_layout(span: float = SPAN, spacing: float = D_SEP) -> Layout:
    k = int(round(span / spacing))
    xs = np.linspace(0, span, k + 1)[1:-1]
    lines = [np.array([[x, 0.0], [x, span]]) for x in xs] + [np.array([[0.0, y], [span, y]]) for y in xs]
    return Layout("B", lines)


def traced_layout(name: str, res: Resultants, tensor: np.ndarray, span: float = SPAN,
                  d_sep: float = D_SEP, degenerate_tol: float = 0.05, magnitude_tol: float = 0.0,
                  min_len_frac: float = 1.5, symmetry: str | None = None) -> Layout:
    xy = res.centroid[:, :2]
    b = (0.0, span, 0.0, span)
    lines = []
    for fam in ("max", "min"):
        lines += trace_lines(xy, tensor, b, d_sep=d_sep, which=fam, degenerate_tol=degenerate_tol,
                             magnitude_tol=magnitude_tol, min_len_frac=min_len_frac,
                             symmetry=symmetry)
    # a line lying along a pinned edge adds nothing for buckling (the edge is already restrained) and
    # only takes rib volume: drop lines whose median distance to the boundary is below EDGE_CLEAR
    # (see docs/decisions.md; author to confirm)
    keep = [ln for ln in lines
            if np.median(np.minimum.reduce([ln[:, 0], span - ln[:, 0], ln[:, 1], span - ln[:, 1]])) >= EDGE_CLEAR]
    return Layout(name, keep)


# Tracing settings (see docs/decisions.md; author to confirm):
#  C - the gravity membrane field is strong everywhere but close to isotropic near the crown, so only
#      a low anisotropy cut-off (0.02) is used.
#  D, Dx - the moment field is concentrated near edges and corners; below 15 % of its maximum the
#      principal directions are numerical noise, so lines stop there (magnitude_tol = 0.15).
#  Symmetry: the shell, gravity and the four-direction envelope are D4-symmetric, so C and D are
#      traced on the group-averaged field with whole orbits (an asymmetric rib set would bias the
#      buckling comparison). The +x case alone is only mirror-symmetric about y = L/2.
TRACE_C = dict(degenerate_tol=0.02, magnitude_tol=0.0, symmetry="D4")
TRACE_D = dict(degenerate_tol=0.05, magnitude_tol=0.15, symmetry="D4")
TRACE_DX = dict(degenerate_tol=0.05, magnitude_tol=0.15, symmetry="mirror_y")


def layout_C(res_g: Resultants, **kw) -> Layout:
    return traced_layout("C", res_g, res_g.n, **{**TRACE_C, **kw})


def layout_D(res: dict[str, Resultants], **kw) -> Layout:
    return traced_layout("D", res["+x"], envelope(res, "m"), **{**TRACE_D, **kw})


def layout_Dx(res: dict[str, Resultants], **kw) -> Layout:
    return traced_layout("Dx", res["+x"], abs_tensor(res["+x"].m), **{**TRACE_DX, **kw})


def crossings(lines: list[np.ndarray]) -> int:
    """Number of points where two different ribs cross (end points on the boundary excluded)."""
    geoms = [LineString(ln) for ln in lines]
    k = 0
    for i in range(len(geoms)):
        for j in range(i + 1, len(geoms)):
            g = geoms[i].intersection(geoms[j])
            if g.is_empty:
                continue
            k += len(getattr(g, "geoms", [g])) if g.geom_type in ("Point", "MultiPoint") else 1
    return k


def ribbed_model(layout: Layout, surface, volume: float, h: float = 0.4,
                 t_base: float = T_BASE_RIBBED, width: float = RIB_WIDTH) -> tuple[Model, dict]:
    """Mesh the layout, set the rib depth for equal ``volume`` and return the model and its metrics."""
    mesh = build_mesh(surface, h, ribs=layout.lines)
    probe = Model(mesh, t_base, ribs=[Rib(i, width, 1.0) for i in range(len(layout.lines))], surface=surface)
    a_s = mesh.area()
    length = sum(probe.rib_length(i) for i in range(len(layout.lines)))
    k = crossings(layout.lines)
    depth = (volume - a_s * t_base) / (width * (length - k * width))
    if depth <= 0:
        raise ValueError(f"volume {volume:.3f} m3 is below the base shell volume {a_s * t_base:.3f} m3")
    ribs = [Rib(i, width, depth) for i in range(len(layout.lines))]
    m = Model(mesh, t_base, ribs=ribs, surface=surface, normals=probe.normals)
    st = line_stats(layout.lines)
    metrics = {"rib_length_m": length, "crossings": k, "rib_depth_mm": depth * 1e3,
               "min_radius_m": st["min_radius"], "n_ribs": len(layout.lines),
               "volume_m3": a_s * t_base + width * depth * (length - k * width)}
    return m, metrics
