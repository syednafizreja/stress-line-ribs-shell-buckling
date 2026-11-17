"""Geometrically nonlinear analysis of the imperfect shell (GNIA, linear elastic material).

1. Linear buckling (LBA) under the load case of interest; *NODE PRINT writes the mode shapes at
   the original (mid-surface) nodes.
2. The mode-1 shape, scaled so its largest translation equals the amplitude w0 (L/300, L/500,
   L/1000 of the span), is added to the node coordinates. Its sign is chosen so the peak moves
   down (the weaker side for a gravity-loaded funicular shell); both signs can be run.
3. *STATIC, NLGEOM with load control: the full load is lambda_max x (G + a_h lateral), applied
   by automatic incrementation. CalculiX 2.23 has no arc-length (Riks) option for this input, so
   the limit load is taken as the last converged load level when the increment size collapses.
   The load-displacement history of the node with the largest mode amplitude is kept, so a
   genuine limit point (flat curve) can be told apart from a solver stop (steep curve).

Knock-down factor alpha = lambda_GNIA / lambda_LBA.
"""
from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from cecore.env import run_ccx

from .constants import A_H, G
from .fe import Model, _grav, _model_cards, _set, lateral_dir

_F = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"


def _loads(m: Model, factor: float, angle_deg: float, a_h: float) -> list[str]:
    grav = lambda g, d: [ln for ln in _grav(g, d) if m.ribs or "RIBS" not in ln]  # noqa: E731
    s = grav(factor * G, (0.0, 0.0, -1.0))
    if a_h:
        s += grav(factor * a_h * G, lateral_dir(angle_deg))[1:]
    return s


def _blocks(dat: Path | str, what: str = "displacements") -> list[tuple[float, dict[int, np.ndarray]]]:
    txt = dat if isinstance(dat, str) else dat.read_text()
    parts = re.split(rf"\s*{what} \(vx,vy,vz\) for set \S+ and time\s+({_F})[^\n]*\n", txt)
    out = []
    row = re.compile(rf"^\s*(\d+)\s+({_F})\s+({_F})\s+({_F})\s*$", re.M)
    for k in range(1, len(parts), 2):
        t = float(parts[k])
        out.append((t, {int(a): np.array([float(b), float(c), float(d)])
                        for a, b, c, d in row.findall(parts[k + 1])}))
    return out


@dataclass
class LBAResult:
    factor: float
    mode: dict[int, np.ndarray]


def lba_mode(m: Model, wd: Path, name: str, angle_deg: float = 0.0, a_h: float = A_H,
             threads: int = 8) -> LBAResult:
    s = _model_cards(m) + ["*STEP", "*BUCKLE", "5"] + _loads(m, 1.0, angle_deg, a_h)
    s += ["*NODE PRINT,NSET=ALLN", "U", "*END STEP", ""]
    wd.mkdir(parents=True, exist_ok=True)
    (wd / f"{name}.inp").write_text("\n".join(s))
    r = run_ccx(name, wd, threads=threads)
    if r.returncode != 0:
        raise RuntimeError((r.stdout + r.stderr)[-800:])
    txt = (wd / f"{name}.dat").read_text()
    blk = txt.split("B U C K L I N G   F A C T O R")[-1]
    facs = [float(x) for x in re.findall(r"^\s+\d+\s+([-+0-9.Ee]+)\s*$", blk.split("displacements")[0], re.M)]
    # the .dat holds the static pre-step block first; the mode shapes follow the factor table
    modes = _blocks(blk)
    if len(modes) != len(facs):
        raise RuntimeError(f"LBA {name}: {len(facs)} factors but {len(modes)} mode blocks")
    k = next(i for i, f in enumerate(facs) if f > 0)
    return LBAResult(facs[k], modes[k][1])


def imperfect(m: Model, mode: dict[int, np.ndarray], w0: float, flip: bool = False) -> Model:
    nodes = m.mesh.nodes
    amp = max(np.linalg.norm(v) for v in mode.values())
    ids = list(mode)
    # sign: the node of largest mode amplitude moves down (inward for this upward-convex shell)
    peak = ids[int(np.argmax([np.linalg.norm(mode[n]) for n in ids]))]
    sign = -1.0 if mode[peak][2] > 0 else 1.0
    if flip:
        sign = -sign
    new = {n: tuple(np.add(p, sign * w0 / amp * mode.get(n, np.zeros(3)))) for n, p in nodes.items()}
    return dataclasses.replace(m, mesh=dataclasses.replace(m.mesh, nodes=new))


@dataclass
class GNIAResult:
    lam_limit: float
    reached_full: bool
    history: list[tuple[float, float]]   # (load factor, |u| at the monitor node)
    monitor: int


def gnia(m: Model, wd: Path, name: str, lam_max: float, angle_deg: float = 0.0, a_h: float = A_H,
         monitor: int | None = None, threads: int = 8, inc0: float = 0.05) -> GNIAResult:
    mon = monitor or next(iter(m.mesh.nodes))
    s = _model_cards(m) + _set("MON", [mon])
    s += ["*STEP,NLGEOM,INC=1000", "*STATIC", f"{inc0},1.0,1e-5,{inc0}"] + _loads(m, lam_max, angle_deg, a_h)
    s += ["*NODE PRINT,NSET=MON", "U", "*END STEP", ""]
    wd.mkdir(parents=True, exist_ok=True)
    (wd / f"{name}.inp").write_text("\n".join(s))
    run_ccx(name, wd, threads=threads)   # a stop at the limit point returns non-zero by design
    hist = [(t * lam_max, float(np.linalg.norm(d[mon]))) for t, d in _blocks(wd / f"{name}.dat") if mon in d]
    if not hist:
        raise RuntimeError(f"GNIA {name}: no converged increment")
    return GNIAResult(hist[-1][0], hist[-1][0] >= lam_max * (1 - 1e-9), hist, mon)
