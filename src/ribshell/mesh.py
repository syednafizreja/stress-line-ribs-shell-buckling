"""Gmsh mesh of the plan square, lifted onto the form-found surface.

Eight-node quadrilaterals (Gmsh ``SecondOrderIncomplete``) map one-to-one onto CalculiX S8R:
corners first, then mid-side nodes in the same order. Rib curves are fragmented with the plan
surface, so the shell mesh has nodes along every rib and at every rib crossing. Each rib's 3-node
line elements become B32R beams. Mid-side nodes are lifted too, so elements follow the curved
surface.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import gmsh
import numpy as np

from .form import Surface

_QUAD8, _LINE3 = 16, 8  # Gmsh element type codes


@dataclass
class ShellMesh:
    nodes: dict[int, tuple[float, float, float]]
    quads: dict[int, tuple[int, ...]]                              # element -> 8 nodes, S8R order
    ribs: dict[int, list[tuple[int, int, int]]] = field(default_factory=dict)  # rib -> (end, mid, end)
    edge_nodes: set[int] = field(default_factory=set)              # plan boundary: pinned

    def centroids(self) -> dict[int, np.ndarray]:
        return {e: np.mean([self.nodes[n] for n in c[:4]], axis=0) for e, c in self.quads.items()}

    def element_areas(self) -> dict[int, float]:
        """Mid-surface area of each element from its 4 corners (two triangles)."""
        out = {}
        for e, c in self.quads.items():
            p = [np.array(self.nodes[n]) for n in c[:4]]
            out[e] = 0.5 * (np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
                            + np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0])))
        return out

    def area(self) -> float:
        return float(sum(self.element_areas().values()))

    def corner_quality(self) -> np.ndarray:
        """Smallest sin(corner angle) of each element in plan (1 = rectangle, <= 0 = degenerate)."""
        out = []
        for c in self.quads.values():
            p = np.array([self.nodes[n][:2] for n in c[:4]])
            a, b = np.roll(p, -1, 0) - p, np.roll(p, 1, 0) - p
            cr = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
            out.append(float(np.min(cr / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)))))
        return np.array(out)


GMSH_OPTS: dict[str, float] = {
    "Mesh.Algorithm": 8,               # frontal-Delaunay for quads
    "Mesh.RecombineAll": 1,
    "Mesh.RecombinationAlgorithm": 2,  # simple full-quad (blossom left 180-deg corners next to ribs)
    "Mesh.Smoothing": 20,
    "Mesh.SubdivisionAlgorithm": 1,    # all-quad: every element is split in four, so the edge is ~h/2
}
#: Quality gate: smallest sine of a corner angle in plan (scaled Jacobian at the corners).
#: The blossom recombination gave corners up to 179.5 deg next to embedded rib curves, which
#: CalculiX rejects (nonpositive Jacobian of the expanded C3D20R); 0.2 is a corner of 11.5-168.5 deg.
MIN_QUALITY = 0.2


def build_mesh(surface: Surface, h: float, ribs: list[np.ndarray] | None = None) -> ShellMesh:
    """Mesh the plan [0, L]^2 with target size ``h`` and lift it onto ``surface``.

    Args:
        surface: Form-found surface.
        h: Target element size in plan [m].
        ribs: Plan polylines (k x 2 arrays). Each is embedded as a spline (a line if k = 2).
    """
    L = surface.span
    ribs = ribs or []
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("plan")
        occ = gmsh.model.occ
        sq = occ.addRectangle(0, 0, 0, L, L)
        curves = []
        for r in ribs:
            pts = [occ.addPoint(float(x), float(y), 0.0) for x, y in r]
            curves.append(occ.addSpline(pts) if len(pts) > 2 else occ.addLine(*pts))
        out_map = None
        if curves:
            _, out_map = occ.fragment([(2, sq)], [(1, c) for c in curves])
        occ.synchronize()
        surfaces = [t for _, t in gmsh.model.getEntities(2)]
        if curves:
            for i in range(len(curves)):
                gmsh.model.addPhysicalGroup(1, [t for d, t in out_map[1 + i] if d == 1], tag=100 + i)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25 * h)
        for k, v in GMSH_OPTS.items():
            gmsh.option.setNumber(k, v)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)   # 8-node quads
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.model.mesh.generate(2)

        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)
        zs = surface.z(xyz[:, 0], xyz[:, 1])
        nodes = {int(t): (float(p[0]), float(p[1]), float(z)) for t, p, z in zip(tags, xyz, zs)}

        quads: dict[int, tuple[int, ...]] = {}
        for s in surfaces:
            types, etags, enodes = gmsh.model.mesh.getElements(2, s)
            for ty, et, en in zip(types, etags, enodes):
                if ty != _QUAD8:
                    raise RuntimeError(f"element type {ty} in shell mesh; expected 8-node quads only")
                for e, conn in zip(et, en.reshape(-1, 8)):
                    c = tuple(int(n) for n in conn)
                    p = np.array([nodes[n][:2] for n in c[:4]])
                    signed = 0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
                    if signed < 0:  # counter-clockwise in plan, so the S8R normal points up (+z)
                        c = (c[0], c[3], c[2], c[1], c[7], c[6], c[5], c[4])
                    quads[int(e)] = c

        rib_elems: dict[int, list[tuple[int, int, int]]] = {}
        for i in range(len(curves)):
            segs = []
            for t in gmsh.model.getEntitiesForPhysicalGroup(1, 100 + i):
                types, _, enodes = gmsh.model.mesh.getElements(1, t)
                for ty, en in zip(types, enodes):
                    if ty == _LINE3:  # Gmsh order: end, end, mid -> B32R order: end, mid, end
                        segs += [(int(a), int(m), int(b)) for a, b, m in en.reshape(-1, 3)]
            rib_elems[i] = segs
    finally:
        gmsh.finalize()
    used = {n for c in quads.values() for n in c}
    nodes = {n: p for n, p in nodes.items() if n in used}
    tol = 1e-6 * L
    edge = {n for n, (x, y, _) in nodes.items() if min(x, y, L - x, L - y) < tol}
    mesh = ShellMesh(nodes, quads, rib_elems, edge)
    q = mesh.corner_quality()
    if q.min() < MIN_QUALITY:
        raise RuntimeError(f"mesh quality gate: {int((q < MIN_QUALITY).sum())} elements with "
                           f"sin(corner) < {MIN_QUALITY} (min {q.min():.3f}), h = {h}")
    return mesh
