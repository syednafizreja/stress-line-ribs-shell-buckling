"""Variants A-D (and Dx) at equal volume V60: buckling factors and fabrication metrics.

For every variant:
  lambda_g        gravity alone (lowest positive buckling factor of 10 modes)
  lambda_c(theta) G + a_h = 0.2 g lateral at theta = 0, 45, ..., 315 deg; lambda_c = min over theta
All eight directions are run for every variant, even where the layout is D4-symmetric and two would
do; the spread across symmetric directions is reported as a check of the mesh and solver symmetry.

Inputs : results/layouts.json (from s01)
Outputs: results/variants_V60.json, results/variants_V60.txt
Usage  : python study/s02_variants.py [--h 0.4]
"""
from __future__ import annotations

import argparse
import json
import tempfile
import time
from pathlib import Path

import numpy as np

from ribshell.constants import T_REF
from ribshell.fe import Model, run_buckle
from ribshell.form import shell_surface
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ANGLES = tuple(range(0, 360, 45))


def lowest_positive(f: list[float]) -> float:
    pos = [x for x in f if x > 0]
    if not pos:
        raise RuntimeError(f"no positive buckling factor in {f}")
    return min(pos)


def evaluate(m: Model, tag: str, wd: Path) -> dict:
    lg = lowest_positive(run_buckle(m, wd, f"{tag}_g", "gravity").factors)
    lc = {a: lowest_positive(run_buckle(m, wd, f"{tag}_c{a}", "combined", a).factors) for a in ANGLES}
    worst = min(lc, key=lc.get)
    return {"lambda_g": lg, "lambda_c": lc[worst], "worst_angle": worst,
            "lambda_c_by_angle": lc, "sway_ratio": lc[worst] / lg}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.4)
    h = ap.parse_args().h
    s = shell_surface()
    lay = json.load(open(RES / "layouts.json"))
    a_mesh = build_mesh(s, h)
    v60 = a_mesh.area() * T_REF
    out: dict[str, dict] = {"h": h, "V60_m3": v60, "variants": {}}
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        t0 = time.perf_counter()
        mA = Model(a_mesh, T_REF, surface=s)
        out["variants"]["A"] = {**evaluate(mA, "A", wd), "volume_m3": mA.volume(), "rib_length_m": 0.0,
                                "n_ribs": 0, "crossings": 0}
        print("A", round(out["variants"]["A"]["lambda_c"], 2), flush=True)
        for name in ("B", "C", "D", "Dx"):
            L = Layout(name, [np.array(ln) for ln in lay[name]])
            m, met = ribbed_model(L, s, v60, h=h)
            met["volume_check_m3"] = m.volume()
            out["variants"][name] = {**evaluate(m, name, wd), **met}
            print(name, round(out["variants"][name]["lambda_c"], 2), flush=True)
        out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(RES / "variants_V60.json", "w"), indent=1)
    rows = [f"Variants at V60 = {v60:.3f} m3, h = {h} m, a_h = 0.2 g (lambda: lowest positive of 10 modes)",
            "var  lam_g    lam_c   worst  lam_c/lam_g  rib_L[m]  depth[mm]  cross  V_model[m3]  sym.spread"]
    for k, r in out["variants"].items():
        lc = r["lambda_c_by_angle"]
        # D4 orbits of the load direction: {0,90,180,270} and {45,135,225,315}
        spread = max(np.ptp([lc[a] for a in (0, 90, 180, 270)]), np.ptp([lc[a] for a in (45, 135, 225, 315)]))
        rows.append(f"{k:4s} {r['lambda_g']:7.2f} {r['lambda_c']:7.2f}  {r['worst_angle']:4d}   {r['sway_ratio']:6.3f}"
                    f"    {r.get('rib_length_m', 0):7.1f}  {r.get('rib_depth_mm', 0):8.0f}  {r['crossings']:5d}"
                    f"   {r.get('volume_check_m3', r['volume_m3']):8.3f}   {spread / r['lambda_c'] * 100:6.2f}%")
    rows.append(f"wall {out['wall_s']:.0f} s")
    (RES / "variants_V60.txt").write_text("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
