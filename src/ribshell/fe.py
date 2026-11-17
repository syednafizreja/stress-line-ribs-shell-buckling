"""CalculiX decks and result readers for the P2 shell.

Conventions, checked on ccx 2.23 (``verify/ccx_capabilities.py``, manual sections 6.2, 7.3, 7.5):

* S8R elements are counter-clockwise in plan, so the shell normal points up (+z).
* Ribs are B32R beams on shell nodes. n2 is the surface normal and n1 = n2 x t, so a
  positive OFFSET2 hangs the rib below the shell (manual 7.3, checked against the .frd).
* Buckling definitions (manual 7.5: earlier loads are dropped unless the buckle step is
  ``*STEP, PERTURBATION``):

  - ``"lateral"``: gravity G is a *STATIC preload and the buckle step is PERTURBATION with a_h.
    lambda_h scales a_h only, so the critical lateral acceleration is a_cr = lambda_h * a_h.
  - ``"combined"``: no preload; lambda_c scales G + a_h together (the scope's literal wording).
  - ``"gravity"``: no preload; lambda_g scales G alone.
"""
from __future__ import annotations

import math
import re
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from cecore.env import run_ccx

from .constants import A_H, E, G, NU, RHO
from .mesh import ShellMesh

ITEMS_PER_LINE = 16


@dataclass
class Rib:
    """One rib: a B32R chain along mesh rib ``rib_id``, rectangular width x depth, hanging below."""

    rib_id: int
    width: float
    depth: float


@dataclass
class Model:
    mesh: ShellMesh
    thickness: float | dict[int, float]      # uniform, or element -> thickness
    ribs: list[Rib] = field(default_factory=list)
    normals: dict[int, np.ndarray] = field(default_factory=dict)  # node -> unit surface normal (ribs)
    surface: object | None = None            # form.Surface; fills rib normals when ``normals`` is empty
    supports: set[int] | None = None         # pinned nodes; default: every plan-boundary node
    e: float = E
    nu: float = NU
    rho: float = RHO

    def __post_init__(self) -> None:
        if self.ribs and not self.normals:
            if self.surface is None:
                raise ValueError("rib models need `surface` (or explicit `normals`) for the rib orientation")
            mids = [mid for r in self.ribs for _, mid, _ in self.mesh.ribs[r.rib_id]]
            xy = np.array([self.mesh.nodes[n][:2] for n in mids])
            nv = self.surface.normal(xy[:, 0], xy[:, 1])
            self.normals = {n: nv[k] for k, n in enumerate(mids)}

    def t_of(self, e: int) -> float:
        return self.thickness if isinstance(self.thickness, float) else self.thickness[e]

    def volume(self) -> float:
        """Shell (mid-surface area x thickness) plus ribs (width x depth x chord length)."""
        areas = self.mesh.element_areas()
        v = sum(a * self.t_of(e) for e, a in areas.items())
        for r in self.ribs:
            v += r.width * r.depth * self.rib_length(r.rib_id)
        return v

    def rib_length(self, rib_id: int) -> float:
        p = self.mesh.nodes
        return float(sum(np.linalg.norm(np.subtract(p[a], p[m])) + np.linalg.norm(np.subtract(p[m], p[b]))
                         for a, m, b in self.mesh.ribs[rib_id]))


def _set(name: str, ids, kind: str = "NSET") -> list[str]:
    ids = sorted(ids)
    return [f"*{kind},{kind}={name}"] + [",".join(map(str, ids[k:k + ITEMS_PER_LINE]))
                                         for k in range(0, len(ids), ITEMS_PER_LINE)]


def _model_cards(m: Model) -> list[str]:
    mesh = m.mesh
    s = ["*NODE"] + [f"{n},{x:.10g},{y:.10g},{z:.10g}" for n, (x, y, z) in mesh.nodes.items()]
    s += ["*ELEMENT,TYPE=S8R,ELSET=SHELL"] + [f"{e}," + ",".join(map(str, c)) for e, c in mesh.quads.items()]
    s += ["*MATERIAL,NAME=C", "*ELASTIC", f"{m.e},{m.nu}", "*DENSITY", f"{m.rho}"]
    # shell sections, one ELSET per distinct thickness (0.01 mm resolution)
    groups: dict[float, list[int]] = defaultdict(list)
    for e in mesh.quads:
        groups[round(m.t_of(e), 5)].append(e)
    for k, (t, els) in enumerate(sorted(groups.items())):
        s += _set(f"T{k}", els, "ELSET") + [f"*SHELL SECTION,ELSET=T{k},MATERIAL=C", f"{t}"]
    # ribs: one ELSET per beam element, because n1 changes along a curved rib on a curved surface
    beam_id = max(mesh.quads) + 1
    rib_sets = []
    # node -> one element containing it (local shell thickness under a rib); built per deck, never cached
    shell_at = {n: e for e, c in mesh.quads.items() for n in c} if m.ribs else {}
    for r in m.ribs:
        for a, mid, b in mesh.ribs[r.rib_id]:
            pa, pb = np.array(mesh.nodes[a]), np.array(mesh.nodes[b])
            t = (pb - pa) / np.linalg.norm(pb - pa)
            n2 = m.normals[mid]
            n1 = np.cross(n2, t)
            n1 /= np.linalg.norm(n1)
            name = f"B{beam_id}"
            t_sh = m.t_of(shell_at[mid])
            off2 = (t_sh / 2 + r.depth / 2) / r.depth  # rib axis below the mid-surface
            s += [f"*ELEMENT,TYPE=B32R,ELSET={name}", f"{beam_id},{a},{mid},{b}",
                  f"*BEAM SECTION,ELSET={name},MATERIAL=C,SECTION=RECT,OFFSET2={off2:.6f}",
                  f"{r.width},{r.depth}", f"{n1[0]:.8f},{n1[1]:.8f},{n1[2]:.8f}"]
            rib_sets.append(name)
            beam_id += 1
    if rib_sets:
        s += ["*ELSET,ELSET=RIBS"] + [",".join(rib_sets[k:k + 8]) for k in range(0, len(rib_sets), 8)]
    sup = mesh.edge_nodes if m.supports is None else m.supports
    s += _set("EDGE", sup) + ["*BOUNDARY", "EDGE,1,3"]
    s += _set("ALLN", mesh.nodes)
    return s


def _grav(g: float, d) -> list[str]:
    els = ["SHELL", "RIBS"]
    return ["*DLOAD"] + [f"{s},GRAV,{g:.10g},{d[0]:.10g},{d[1]:.10g},{d[2]:.10g}" for s in els]


def lateral_dir(angle_deg: float) -> tuple[float, float, float]:
    a = math.radians(angle_deg)
    return (math.cos(a), math.sin(a), 0.0)


def buckle_deck(m: Model, kind: str, angle_deg: float = 0.0, a_h: float = A_H, modes: int = 10,
                solver: str | None = None, frd: bool = False) -> str:
    """Deck for one buckling analysis. ``kind`` is "lateral", "combined" or "gravity"."""
    s = _model_cards(m)
    # GRAV on the RIBS set only when ribs exist (an undefined set is a ccx error)
    grav = lambda g, d: [ln for ln in _grav(g, d) if m.ribs or "RIBS" not in ln]  # noqa: E731
    buck = f"*BUCKLE,SOLVER={solver}" if solver else "*BUCKLE"
    down, lat = (0.0, 0.0, -1.0), lateral_dir(angle_deg)
    if kind == "lateral":
        s += ["*STEP", "*STATIC"] + grav(G, down) + ["*END STEP"]
        s += ["*STEP,PERTURBATION", buck, f"{modes}"] + grav(a_h * G, lat)
    elif kind == "combined":
        s += ["*STEP", buck, f"{modes}"] + grav(G, down) + grav(a_h * G, lat)
    elif kind == "gravity":
        s += ["*STEP", buck, f"{modes}"] + grav(G, down)
    else:
        raise ValueError(kind)
    s += ["*NODE FILE", "U", "*END STEP", ""]
    return "\n".join(s)


def read_factors(dat: Path) -> list[float]:
    txt = dat.read_text()
    blk = txt.split("B U C K L I N G   F A C T O R")[-1] if "B U C K L I N G" in txt else ""
    return [float(x) for x in re.findall(r"^\s+\d+\s+([-+0-9.Ee]+)\s*$", blk, re.M)]


@dataclass
class BuckleResult:
    factors: list[float]
    wall_s: float
    solver: str
    workdir: Path


def run_buckle(m: Model, workdir: Path, name: str, kind: str, angle_deg: float = 0.0,
               a_h: float = A_H, modes: int = 10, threads: int = 4, frd: bool = False,
               retries: int = 1) -> BuckleResult:
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / f"{name}.inp").write_text(buckle_deck(m, kind, angle_deg, a_h, modes, frd=frd))
    t0 = time.perf_counter()
    r = run_ccx(name, workdir, threads=threads)
    for _ in range(retries):   # one retry for transient OS failures (out-of-memory status, file locks)
        if r.returncode == 0:
            break
        time.sleep(5)
        r = run_ccx(name, workdir, threads=threads)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        raise RuntimeError(f"ccx failed ({r.returncode}): {(r.stdout + r.stderr)[-600:]}")
    log = r.stdout + r.stderr
    solver = next((s for s in ("PARDISO", "PaStiX", "SPOOLES", "TAUCS") if s.lower() in log.lower()), "?")
    res = BuckleResult(read_factors(workdir / f"{name}.dat"), wall, solver, workdir)
    if not frd:   # keep the work folder small during long studies
        for ext in (".inp", ".dat", ".sta", ".cvg", ".12d", ".frd"):
            (workdir / f"{name}{ext}").unlink(missing_ok=True)
    return res
