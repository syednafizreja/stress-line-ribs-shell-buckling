"""Linear static load cases and shell stress resultants (n, m) per element.

CalculiX expands each S8R into a C3D20R through the thickness and integrates it with 2 x 2 x 2
Gauss points. ``*EL PRINT, S`` writes the 8 integration-point stresses in global axes. Under
Kirchhoff-Love/Mindlin kinematics the in-plane stress is linear through the thickness, so with the
two Gauss layers at zeta = +-1/sqrt(3) (z = +-t / (2 sqrt 3)):

    n_ab = t * sigma_mid,           sigma_mid = (sigma_top + sigma_bot) / 2
    m_ab = sqrt(3) t^2 (sigma_top - sigma_bot) / 12

both exact. Sign convention: m > 0 puts the +normal (top, upper) face in tension. Which four
integration points form the top layer is fixed by ``verify/check_fields.py`` (a sagging plate).

The tensors are projected on a local tangent frame e1 = unit(1, 0, z_x), e2 = n x e1, which follows
the plan x axis; for this height-field shell the tracer works in plan with these components.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from cecore.env import run_ccx

from .constants import A_H, G
from .fe import Model, _grav, _model_cards, lateral_dir

TOP_IPS = (4, 5, 6, 7)   # 0-based: integration points 5-8 lie on the +normal side (verified)
BOT_IPS = (0, 1, 2, 3)


@dataclass
class Resultants:
    """Element-centroid stress resultants in the local tangent frame (N/m and N m/m)."""

    centroid: np.ndarray   # (ne, 3)
    n: np.ndarray          # (ne, 2, 2)
    m: np.ndarray          # (ne, 2, 2)
    elems: np.ndarray      # (ne,) element ids


def static_cases(a_h: float = A_H) -> dict[str, list[tuple[float, tuple[float, float, float]]]]:
    """Gravity alone and G + a_h in +x, -x, +y, -y (the envelope directions)."""
    down = (0.0, 0.0, -1.0)
    cases = {"G": [(G, down)]}
    for name, ang in (("+x", 0), ("+y", 90), ("-x", 180), ("-y", 270)):
        cases[name] = [(G, down), (a_h * G, lateral_dir(ang))]
    return cases


def static_deck(m: Model, cases) -> str:
    s = _model_cards(m)
    for loads in cases.values():
        s += ["*STEP", "*STATIC"]
        for k, (g, d) in enumerate(loads):
            lines = [ln for ln in _grav(g, d) if m.ribs or "RIBS" not in ln]
            if k == 0:
                lines[0] = "*DLOAD,OP=NEW"
            s += lines
        s += ["*NODE FILE", "U", "*EL PRINT,ELSET=SHELL", "S", "*END STEP"]
    return "\n".join(s) + "\n"


_FLOAT = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"


def read_ip_stresses(dat: Path, n_steps: int) -> list[dict[int, np.ndarray]]:
    """Per step: element -> (8, 6) array of sxx, syy, szz, sxy, sxz, syz at the integration points."""
    txt = dat.read_text()
    blocks = re.split(r"\s*stresses \(elem, integ\.pnt\.,sxx,syy,szz,sxy,sxz,syz\)[^\n]*\n", txt)[1:]
    if len(blocks) != n_steps:
        raise RuntimeError(f"expected {n_steps} stress blocks in {dat.name}, found {len(blocks)}")
    out = []
    row = re.compile(rf"^\s*(\d+)\s+(\d+)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})(?:\s+\S+)?[ \t]*$", re.M)
    for b in blocks:
        acc: dict[int, list] = {}
        for mt in row.finditer(b):
            acc.setdefault(int(mt.group(1)), []).append([float(mt.group(k)) for k in range(3, 9)])
        out.append({e: np.array(v) for e, v in acc.items()})
    return out


def _voigt_to_tensor(s: np.ndarray) -> np.ndarray:
    """(..., 6) sxx syy szz sxy sxz syz -> (..., 3, 3)."""
    t = np.empty(s.shape[:-1] + (3, 3))
    t[..., 0, 0], t[..., 1, 1], t[..., 2, 2] = s[..., 0], s[..., 1], s[..., 2]
    t[..., 0, 1] = t[..., 1, 0] = s[..., 3]
    t[..., 0, 2] = t[..., 2, 0] = s[..., 4]
    t[..., 1, 2] = t[..., 2, 1] = s[..., 5]
    return t


def tangent_frame(surface, x, y) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    zx, _ = surface.grad(x, y)
    e1 = np.stack([np.ones_like(zx), np.zeros_like(zx), zx], axis=-1)
    e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    nrm = surface.normal(x, y)
    e2 = np.cross(nrm, e1)
    return e1, e2, nrm


def resultants(m: Model, ip: dict[int, np.ndarray], surface) -> Resultants:
    elems = np.array(sorted(ip))
    cen = m.mesh.centroids()
    c = np.array([cen[e] for e in elems])
    t = np.array([m.t_of(int(e)) for e in elems])
    s = np.stack([ip[e] for e in elems])                    # (ne, 8, 6)
    top = _voigt_to_tensor(s[:, TOP_IPS].mean(axis=1))
    bot = _voigt_to_tensor(s[:, BOT_IPS].mean(axis=1))
    e1, e2, _ = tangent_frame(surface, c[:, 0], c[:, 1])
    basis = np.stack([e1, e2], axis=1)                       # (ne, 2, 3)
    proj = lambda T: np.einsum("eai,eij,ebj->eab", basis, T, basis)  # noqa: E731
    mid, diff = proj(0.5 * (top + bot)), proj(top - bot)
    n = t[:, None, None] * mid
    mm = np.sqrt(3.0) * t[:, None, None] ** 2 * diff / 12.0
    return Resultants(c, n, mm, elems)


def run_static(m: Model, surface, workdir: Path, name: str, cases=None, threads: int = 4
               ) -> dict[str, Resultants]:
    cases = cases or static_cases()
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / f"{name}.inp").write_text(static_deck(m, cases))
    r = run_ccx(name, workdir, threads=threads)
    if r.returncode != 0:
        raise RuntimeError(f"ccx failed: {(r.stdout + r.stderr)[-600:]}")
    steps = read_ip_stresses(workdir / f"{name}.dat", len(cases))
    return {k: resultants(m, ip, surface) for k, ip in zip(cases, steps)}


def abs_tensor(T: np.ndarray) -> np.ndarray:
    """sum |lambda_i| v_i v_i^T for (..., 2, 2) symmetric tensors."""
    w, v = np.linalg.eigh(T)
    return np.einsum("...ik,...k,...jk->...ij", v, np.abs(w), v)


def envelope(res: dict[str, Resultants], which: str = "m", keys=("+x", "-x", "+y", "-y")) -> np.ndarray:
    """Envelope tensor: sum over the lateral cases of the absolute-principal tensor (scope 5.3)."""
    return sum(abs_tensor(getattr(res[k], which)) for k in keys)
