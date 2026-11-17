"""Robustness of the traced layouts to the mesh of the source field.

C, D and Dx are traced from the linear fields of variant A computed at h = 0.6, 0.4, 0.3, 0.2 m.
For each layout: number of lines, total plan length, and the distance to the h = 0.4 layout
(the one used in the study): mean and maximum over the points of one layout of the distance to
the nearest point of the other, taken both ways (symmetric; max = Hausdorff distance).
The gravity membrane field (C) is large and smooth; the seismic moment field (D, Dx) is small
(|m| ~ 0.03 kN m/m on a funicular form), so its lines may change topology with the mesh.

Usage: python verify/layout_robustness.py   (writes verify/out/layout_robustness.{txt,json})
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ribshell.constants import T_REF
from ribshell.fe import Model
from ribshell.fields import run_static
from ribshell.form import shell_surface
from ribshell.layouts import layout_C, layout_D, layout_Dx
from ribshell.mesh import build_mesh

OUT = Path(__file__).parent / "out"
HS = (0.6, 0.4, 0.3, 0.2)
REF_H = 0.4


def resample(lines, ds=0.05):
    pts = []
    for ln in lines:
        seg = np.linalg.norm(np.diff(ln, axis=0), axis=1)
        s = np.concatenate([[0], np.cumsum(seg)])
        t = np.arange(0, s[-1] + 1e-12, ds)
        pts.append(np.column_stack([np.interp(t, s, ln[:, 0]), np.interp(t, s, ln[:, 1])]))
    return np.vstack(pts) if pts else np.zeros((0, 2))


def distance(a, b):
    pa, pb = resample(a), resample(b)
    if not len(pa) or not len(pb):
        return float("nan"), float("nan")
    da, db = cKDTree(pb).query(pa)[0], cKDTree(pa).query(pb)[0]
    return float(np.concatenate([da, db]).mean()), float(max(da.max(), db.max()))


def main() -> None:
    s = shell_surface()
    lays = {}
    with tempfile.TemporaryDirectory() as td:
        for h in HS:
            res = run_static(Model(build_mesh(s, h), T_REF), s, Path(td), f"A_{int(h*100)}")
            lays[h] = {"C": layout_C(res["G"]).lines, "D": layout_D(res).lines, "Dx": layout_Dx(res).lines}
            print(h, {k: len(v) for k, v in lays[h].items()}, flush=True)
    rows, data = ["layout robustness to the field mesh (reference h = 0.4 m)",
                  "lay  h     n_lines  length[m]  mean_d[m]  hausdorff[m]"], {}
    for k in ("C", "D", "Dx"):
        for h in HS:
            ln = lays[h][k]
            length = sum(np.linalg.norm(np.diff(x, axis=0), axis=1).sum() for x in ln)
            md, hd = distance(ln, lays[REF_H][k])
            data.setdefault(k, {})[str(h)] = {"n_lines": len(ln), "length_m": length, "mean_d_m": md,
                                              "hausdorff_m": hd, "lines": [x.round(4).tolist() for x in ln]}
            rows.append(f"{k:3s}  {h:<4}  {len(ln):7d}  {length:9.1f}  {md:9.3f}  {hd:12.3f}")
    OUT.mkdir(exist_ok=True)
    (OUT / "layout_robustness.txt").write_text("\n".join(rows) + "\n")
    json.dump(data, open(OUT / "layout_robustness.json", "w"))
    print("\n".join(rows))


if __name__ == "__main__":
    main()
