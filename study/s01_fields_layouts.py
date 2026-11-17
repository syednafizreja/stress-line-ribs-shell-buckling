"""Step 1: linear fields on variant A (uniform 60 mm) and the rib layouts B, C, D, Dx.

Writes results/fields_A.npz, results/layouts.json and figures/layouts.png.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ribshell.constants import SPAN, T_REF
from ribshell.fe import Model
from ribshell.fields import envelope, run_static
from ribshell.form import shell_surface
from ribshell.layouts import grid_layout, layout_C, layout_D, layout_Dx
from ribshell.mesh import build_mesh

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"


def main() -> None:
    RES.mkdir(exist_ok=True), FIG.mkdir(exist_ok=True)
    s = shell_surface()
    m = Model(build_mesh(s, 0.4), T_REF)
    with tempfile.TemporaryDirectory() as td:
        res = run_static(m, s, Path(td), "A_static")
    np.savez_compressed(RES / "fields_A.npz", centroid=res["G"].centroid,
                        **{f"n_{k}": r.n for k, r in res.items()}, **{f"m_{k}": r.m for k, r in res.items()})
    layouts = [grid_layout(), layout_C(res["G"]), layout_D(res), layout_Dx(res)]
    json.dump({L.name: [ln.round(4).tolist() for ln in L.lines] for L in layouts},
              open(RES / "layouts.json", "w"))
    # summary numbers for the record
    g = res["G"]
    mg = np.abs(np.linalg.eigvalsh(g.m)).max(axis=1)
    ng = np.abs(np.linalg.eigvalsh(g.n)).max(axis=1)
    mx = np.abs(np.linalg.eigvalsh(res["+x"].m)).max(axis=1)
    print(f"gravity: max|n| {ng.max()/1e3:.1f} kN/m, max|m| {mg.max():.1f} N m/m, "
          f"median |m|/(|n| t) {np.median(mg/(ng*T_REF)):.4f}")
    print(f"G+0.2g x: max|m| {mx.max():.1f} N m/m  ({mx.max()/mg.max():.1f} x gravity)")
    for L in layouts:
        print(f"{L.name}: {len(L.lines)} lines, plan length {sum(np.linalg.norm(np.diff(l,axis=0),axis=1).sum() for l in L.lines):.1f} m")

    # figure: one panel per layout over the field it was traced from
    env = envelope(res, "m")
    bg = {"B": None, "C": np.abs(np.linalg.eigvalsh(g.n)).max(axis=1) / 1e3,
          "D": np.linalg.eigvalsh(env).max(axis=1) / 1e3, "Dx": mx / 1e3}
    lab = {"B": "Orthogonal grid, 1.5 m", "C": "Gravity membrane-force lines",
           "D": "Seismic-envelope moment lines", "Dx": "Moment lines, +x case only"}
    unit = {"C": "max |n| [kN/m]", "D": "envelope max [kN m/m]", "Dx": "max |m| [kN m/m]"}
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans"})
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.4), constrained_layout=True)
    c = g.centroid
    for ax, L in zip(axs, layouts):
        if bg[L.name] is not None:
            tc = ax.tricontourf(c[:, 0], c[:, 1], bg[L.name], levels=14, cmap="Greys")
            fig.colorbar(tc, ax=ax, shrink=0.75, label=unit[L.name])
        for ln in L.lines:
            ax.plot(ln[:, 0], ln[:, 1], color="#c0392b", lw=1.1)
        ax.set_aspect("equal"), ax.set_xlim(0, SPAN), ax.set_ylim(0, SPAN)
        ax.set_title(f"{L.name}: {lab[L.name]}", fontsize=8.5)
        ax.set_xlabel("x [m]")
    axs[0].set_ylabel("y [m]")
    fig.savefig(FIG / "layouts.png", dpi=220)


if __name__ == "__main__":
    main()
