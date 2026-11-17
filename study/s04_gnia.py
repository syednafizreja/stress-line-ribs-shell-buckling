"""Imperfection sensitivity: geometrically nonlinear analysis of the imperfect shell (GNIA).

For each variant (A uniform, C and D ribbed, E zonal optimum) at equal volume V60:
  1. LBA in the governing lateral direction (the worst of 0 and 45 deg; by D4 symmetry these two
     bound all eight for A, C, D and the d4 E design);
  2. mode-1 imperfection with amplitude w0 = L/300, L/500, L/1000 (L = 12 m plan span; the L/500
     value is Claude's reference choice, author to confirm against the code basis used);
     both signs of the mode (the weaker one governs);
  3. NLGEOM load control up to 1.2 lambda_LBA x (G + 0.2 g lateral); the limit is the last
     converged increment (CalculiX has no arc-length option for this input).
Knock-down alpha = lambda_GNIA / lambda_LBA. A flat end of the load-displacement curve (secant
stiffness of the last 10 % of load below 10 % of the initial stiffness) marks a genuine limit point;
otherwise the stop is flagged as "solver stop".

Usage: python study/s04_gnia.py [--h 0.6]   -> results/gnia_h<h>.json, results/gnia_h<h>.txt
"""
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from ribshell.constants import SPAN, T_REF
from ribshell.fe import Model
from ribshell.form import shell_surface
from ribshell.gnia import gnia, imperfect, lba_mode
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
W0 = {"L/300": SPAN / 300, "L/500": SPAN / 500, "L/1000": SPAN / 1000}
WORKERS, THREADS = 3, 4


def models(h: float) -> dict[str, Model]:
    s = shell_surface()
    mesh = build_mesh(s, h)
    v60 = mesh.area() * T_REF
    out = {"A": Model(mesh, T_REF, surface=s)}
    lay = json.load(open(RES / "layouts.json"))
    for n in ("C", "D"):
        out[n], _ = ribbed_model(Layout(n, [np.array(x) for x in lay[n]]), s, v60, h=h)
    eb = RES / "E_best_d4_h0.6.json"
    if eb.exists():
        e = json.load(open(eb))
        des = ZonalDesign(mesh, e["nz"], v60)
        tz = des.zone_thickness(np.asarray(e["zone_thickness_m"]))   # re-projected on this mesh
        out["E"] = Model(mesh, des.element_thickness(tz), surface=s)
    return out


def limit_kind(hist: list[tuple[float, float]]) -> str:
    lam = np.array([p[0] for p in hist])
    u = np.array([p[1] for p in hist])
    if len(lam) < 4 or u[0] <= 0:
        return "undetermined"
    k0 = lam[0] / u[0]
    sel = lam >= 0.9 * lam[-1]
    i = int(np.argmax(sel))
    i = min(i, len(lam) - 2)
    k_end = (lam[-1] - lam[i]) / max(u[-1] - u[i], 1e-15)
    return "limit point" if k_end < 0.1 * k0 else "solver stop"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.6)
    a = ap.parse_args()
    ms = models(a.h)
    wd = Path(tempfile.mkdtemp(prefix="gnia_"))
    lba = {}
    for n, m in ms.items():
        r = {ang: lba_mode(m, wd, f"{n}_lba{ang}", ang, threads=8) for ang in (0, 45)}
        ang = min(r, key=lambda k: r[k].factor)
        lba[n] = (ang, r[ang])
        print(f"{n}: LBA {r[0].factor:.2f} (0)  {r[45].factor:.2f} (45) -> governing {ang}", flush=True)

    jobs = [(n, wk, flip) for n in ms for wk in W0 for flip in (False, True)]

    def run(job):
        n, wk, flip = job
        ang, L = lba[n]
        peak = max(L.mode, key=lambda k: float(L.mode[k] @ L.mode[k]))
        mi = imperfect(ms[n], L.mode, W0[wk], flip=flip)
        tag = f"{n}_{wk.replace('/', '')}_{'m' if flip else 'p'}"
        g = gnia(mi, wd, tag, 1.2 * L.factor, ang, monitor=peak, threads=THREADS)
        return {"variant": n, "w0": wk, "w0_m": W0[wk], "sign": "flipped" if flip else "peak-down",
                "angle_deg": ang, "lambda_LBA": L.factor, "lambda_GNIA": g.lam_limit,
                "alpha": g.lam_limit / L.factor, "reached_full": g.reached_full,
                "kind": limit_kind(g.history), "u_peak_end_m": g.history[-1][1],
                "history": g.history}

    with ThreadPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(run, jobs))
    shutil.rmtree(wd, ignore_errors=True)
    json.dump({"h": a.h, "rows": rows}, open(RES / f"gnia_h{a.h:g}.json", "w"), indent=1)
    lines = [f"GNIA, h = {a.h} m, load control to 1.2 lambda_LBA, linear elastic C30/37",
             "var  w0      sign       ang  lam_LBA  lam_GNIA  alpha   kind"]
    for r in rows:
        lines.append(f"{r['variant']:3s}  {r['w0']:6s}  {r['sign']:9s}  {r['angle_deg']:3d}  {r['lambda_LBA']:7.2f}"
                     f"  {r['lambda_GNIA']:8.2f}  {r['alpha']:.3f}  {r['kind']}")
    lines.append("")
    lines.append("governing (min over sign):  var  w0  lam_GNIA  alpha")
    for n in ms:
        for wk in W0:
            g = min((r for r in rows if r["variant"] == n and r["w0"] == wk), key=lambda r: r["lambda_GNIA"])
            lines.append(f"  {n:3s} {wk:6s} {g['lambda_GNIA']:8.2f}  {g['alpha']:.3f}")
    (RES / f"gnia_h{a.h:g}.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
