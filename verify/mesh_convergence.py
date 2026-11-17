"""Mesh convergence of the buckling factors (brief section 7: three meshes).

Variants A (uniform 60 mm), C and D (ribbed, equal volume) and E (best zonal design, if present)
at h = 0.6, 0.4, 0.3, 0.2 m: lambda_g and lambda_c at 0 and 45 deg. For the three finest meshes the
observed order p and the Richardson estimate lambda_inf are reported; the change from h = 0.3
to 0.2 is the discretisation error quoted in the report.

Usage: python verify/mesh_convergence.py   (writes verify/out/mesh_convergence.{txt,json})
"""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import numpy as np

from ribshell.constants import T_REF
from ribshell.fe import Model, run_buckle
from ribshell.form import shell_surface
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).parent / "out"
HS = (0.6, 0.4, 0.3, 0.2)


def low(f):
    return min(x for x in f if x > 0)


def build(name, s, h, lay, e_best):
    mesh_a = build_mesh(s, h)
    v60 = mesh_a.area() * T_REF   # V60 on this mesh (the area changes by < 0.1 % between meshes)
    if name == "A":
        return Model(mesh_a, T_REF, surface=s)
    if name == "E":
        des = ZonalDesign(mesh_a, e_best["nz"], v60)
        tz = des.zone_thickness(np.asarray(e_best["zone_thickness_m"]))
        return Model(mesh_a, des.element_thickness(tz), surface=s)
    m, _ = ribbed_model(Layout(name, [np.array(x) for x in lay[name]]), s, v60, h=h)
    return m


def richardson(h, f):
    """Observed order and extrapolated value from three meshes (ratio may be non-constant)."""
    (h1, h2, h3), (f1, f2, f3) = h, f
    try:
        r = (f1 - f2) / (f2 - f3)
        p = math.log(abs(r)) / math.log(h1 / h2) if r > 0 else float("nan")
        finf = f3 + (f3 - f2) / ((h2 / h3) ** p - 1) if p == p else float("nan")
    except (ZeroDivisionError, ValueError):
        p, finf = float("nan"), float("nan")
    return p, finf


def main() -> None:
    s = shell_surface()
    lay = json.load(open(ROOT / "results" / "layouts.json"))
    eb = ROOT / "results" / "E_best_d4_h0.6.json"
    e_best = json.load(open(eb)) if eb.exists() else None
    names = ["A", "C", "D"] + (["E"] if e_best else [])
    data = {}
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        for n in names:
            data[n] = {}
            for h in HS:
                m = build(n, s, h, lay, e_best)
                tag = f"{n}_{int(h*100)}"
                r = {"lambda_g": low(run_buckle(m, wd, tag + "g", "gravity", threads=8).factors),
                     "lambda_c0": low(run_buckle(m, wd, tag + "c0", "combined", 0, threads=8).factors),
                     "lambda_c45": low(run_buckle(m, wd, tag + "c45", "combined", 45, threads=8).factors),
                     "n_shell": len(m.mesh.quads)}
                data[n][h] = r
                print(n, h, {k: round(v, 2) for k, v in r.items()}, flush=True)
    lines = ["Mesh convergence of buckling factors (lowest positive of 10 modes)",
             "var key         " + "  ".join(f"h={h:<5}" for h in HS) + "  d(0.3->0.2)  p_obs  lam_inf"]
    for n, d in data.items():
        for k in ("lambda_g", "lambda_c0", "lambda_c45"):
            vals = [d[h][k] for h in HS]
            p, finf = richardson(HS[1:], vals[1:])
            lines.append(f"{n:3s} {k:10s} " + "  ".join(f"{v:7.2f}" for v in vals)
                         + f"   {(vals[-1]-vals[-2])/vals[-1]*100:+6.2f}%   {p:5.2f}  {finf:7.2f}")
    OUT.mkdir(exist_ok=True)
    (OUT / "mesh_convergence.txt").write_text("\n".join(lines) + "\n")
    json.dump({n: {str(h): v for h, v in d.items()} for n, d in data.items()},
              open(OUT / "mesh_convergence.json", "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
