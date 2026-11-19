"""All report figures from the result files (no analysis is run here).

Reads   results/sweep_h<h>.json, variants_h<h>_r1.json, E_best_{d4,free}_h<h>.json, gnia_h<h>.json,
        layouts.json and modes_h<h>.npz (written by s07_exports.py).
Writes  figures/fig1_lambda_vs_volume.png   lambda_c vs concrete volume, tCO2e on the top axis,
                                            marker area = total rib length (fabrication effort)
        figures/fig2_recovery.png           recovery R_X = (lambda_X - lambda_A) / (lambda_E - lambda_A)
        figures/fig3_E_thickness.png        zone thickness of E (D4 and free optimum, V = V60)
        figures/fig4_gnia.png               load-displacement paths with the first-mode imperfection
        figures/fig5_modes.png              plan view of mode 1 under sway for A, C, D and E
        figures/fig6_by_angle.png           lambda_c by sway direction at V = V60
        figures/table_*.csv                 the data behind each figure (table view)
A figure whose input file is missing is skipped with a message.
Usage   python study/s06_figures.py [--h 0.4]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"

# Categorical order fixed per variant (checked with the dataviz palette validator: CVD and
# normal-vision separation pass; three hues are below 3:1 contrast, so every line is also
# direct-labelled and has its own marker).
COL = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a", "D": "#eda100", "Dd": "#e87ba4",
       "Dx": "#008300", "E": "#6250d6"}
MRK = {"A": "o", "B": "s", "C": "^", "D": "D", "Dd": "v", "Dx": "P", "E": "*"}
NAME = {"A": "A uniform", "B": "B grid", "C": "C gravity lines", "D": "D sway-envelope lines",
        "Dd": "Dd diagonal lines", "Dx": "Dx +x sway lines", "E": "E optimised thickness"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.dpi": 220,
    "legend.frameon": False, "lines.linewidth": 2})


def load(name: str):
    p = RES / name
    if not p.exists():
        print(f"skip: {p.name} not found")
        return None
    return json.load(open(p))


def table(name: str, header: list[str], rows: list[list]) -> None:
    with open(FIG / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def fig_sweep(sw: dict) -> None:
    rows = sw["rows"]
    names = [k for k in COL if k in rows[0]["lambda_c"]]
    V = np.array([r["V_m3"] for r in rows])
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ends: list[tuple[float, str]] = []
    for k in names:
        lam = np.array([r["lambda_c"][k] for r in rows])
        ax.plot(V, lam, color=COL[k], zorder=2)
        if k in rows[0].get("fabrication", {}):
            L = np.array([r["fabrication"][k]["rib_length_m"] for r in rows])
            ax.scatter(V, lam, s=12 + L * 0.9, marker=MRK[k], color=COL[k], edgecolor="#fcfcfb",
                       linewidth=1.2, zorder=3)
        else:
            ax.scatter(V, lam, s=40 if k != "E" else 90, marker=MRK[k], color=COL[k],
                       edgecolor="#fcfcfb", linewidth=1.2, zorder=3)
        ends.append((float(lam[-1]), k))
    # Keep the end labels off each other. D and Dd finish within one line height.
    ends.sort()
    gap, shown = 22.0, []
    for y, k in ends:
        yy = y if not shown or y >= shown[-1] + gap else shown[-1] + gap
        shown.append(yy)
        if yy == y:
            ax.annotate(k, (V[-1], y), xytext=(8, 0), textcoords="offset points",
                        va="center", color=INK, fontsize=9, fontweight="bold")
        else:
            ax.annotate(k, (V[-1], y), xytext=(V[-1] + 0.35, yy), textcoords="data",
                        va="center", color=INK, fontsize=9, fontweight="bold",
                        arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.6})
    if "E_free" in rows[0]:
        ax.scatter([V[0]], [rows[0]["E_free"]], marker="*", s=90, facecolor="none", edgecolor=COL["E"],
                   linewidth=1.4, zorder=4, label="E, free 25-zone optimum (V60 only)")
    ax.set_xlabel("Concrete volume V [m³]  (A1–A3 carbon on the top axis)")
    ax.set_ylabel("λ_c, lowest factor over 8 sway directions")
    ecf, rho = sw["ecf_kgco2e_per_kg"]["default"], 2500.0
    top = ax.secondary_xaxis("top", functions=(lambda v: v * rho * ecf / 1e3, lambda c: c * 1e3 / (rho * ecf)))
    top.set_xlabel(f"Embodied carbon A1–A3 [tCO₂e]  ({ecf} kgCO₂e/kg)", color=INK2)
    for L in (50, 100, 150):
        ax.scatter([], [], s=12 + L * 0.9, color="#b9b8b2", label=f"rib length {L} m")
    ax.legend(loc="upper left", fontsize=8, handletextpad=0.4)
    ax.set_xlim(V[0] - 0.3, V[-1] + 0.9)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_lambda_vs_volume.png")
    plt.close(fig)
    table("table_sweep.csv", ["V_over_V60", "V_m3", "tCO2e", *[f"lambda_c_{k}" for k in names]],
          [[r["vratio"], round(r["V_m3"], 3), round(r["tCO2e"]["default"], 2),
            *[round(r["lambda_c"][k], 2) for k in names]] for r in rows])

    rib = [k for k in names if k in rows[0]["recovery"]]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.array([r["vratio"] for r in rows])
    w = 0.8 / len(rib) * 0.1
    for i, k in enumerate(rib):
        R = np.array([r["recovery"][k] for r in rows]) * 100
        ax.bar(x + (i - (len(rib) - 1) / 2) * w, R, width=w * 0.9, color=COL[k], label=NAME[k], zorder=2)
    ax.axhline(70, color=INK, lw=1, ls="--")
    ax.text(x[-1], 90, "hypothesis: D ≥ 70 %", fontsize=8, color=INK, ha="right")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x, [f"{v:.1f}" for v in x])
    ax.set_xlabel("V / V60")
    ax.set_ylabel("Recovery of the E gain over A [%]")
    ax.legend(ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    fig.tight_layout()
    fig.savefig(FIG / "fig2_recovery.png", bbox_inches="tight")
    plt.close(fig)
    table("table_recovery.csv", ["V_over_V60", *[f"R_{k}" for k in rib]],
          [[r["vratio"], *[round(r["recovery"][k], 4) for k in rib]] for r in rows])


def fig_E(h: float) -> None:
    eds = [(s, load(f"E_best_{s}_h{h:g}.json")) for s in ("d4", "free")]
    eds = [(s, e) for s, e in eds if e]
    if not eds:
        return
    fig, axs = plt.subplots(1, len(eds), figsize=(3.6 * len(eds) + 0.8, 3.4))
    axs = np.atleast_1d(axs)
    tz_all = np.concatenate([np.asarray(e["zone_thickness_m"]) for _, e in eds]) * 1e3
    vmin, vmax = tz_all.min(), tz_all.max()
    for ax, (s, e) in zip(axs, eds):
        nz = e["nz"]
        # zonal.zone_of_elements numbers zones ix * nz + iy; transpose to rows = y for imshow
        tz = np.asarray(e["zone_thickness_m"]).reshape(nz, nz).T * 1e3
        im = ax.imshow(tz, origin="lower", cmap="Blues", vmin=vmin, vmax=vmax, extent=(0, 12, 0, 12))
        for iy in range(nz):
            for ix in range(nz):
                c = (ix + 0.5) * 12 / nz, (iy + 0.5) * 12 / nz
                ax.text(*c, f"{tz[iy, ix]:.0f}", ha="center", va="center", fontsize=7.5,
                        color="#fcfcfb" if tz[iy, ix] > (vmin + vmax) / 2 else INK)
        ax.set_title(f"E {'D4-symmetric (6 vars)' if s == 'd4' else 'free (25 vars)'}\nλ_c = {e['objective']:.1f}",
                     fontsize=9)
        ax.set_xlabel("x [m]")
        ax.grid(False)
    axs[0].set_ylabel("y [m]")
    fig.colorbar(im, ax=axs.tolist(), label="shell thickness [mm]", shrink=0.85)
    fig.savefig(FIG / "fig3_E_thickness.png", bbox_inches="tight")
    plt.close(fig)


def fig_gnia(g: dict) -> None:
    rows = g["rows"]
    names = [k for k in COL if any(r["variant"] == k for r in rows)]
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.8))
    for r in rows:
        if r["w0"] != "L/500":
            continue
        hst = np.array(r["history"])
        ls = "-" if r["sign"] == "peak-down" else ":"
        axs[0].plot(np.abs(hst[:, 1]) * 1e3, hst[:, 0] / r["lambda_LBA"], ls=ls, color=COL[r["variant"]],
                    label=f"{r['variant']} ({r['sign']})")
    axs[0].set_xlabel("|u| at the mode-1 peak node [mm]")
    axs[0].set_ylabel("λ / λ_LBA")
    axs[0].set_title("w₀ = L/500, both imperfection signs", fontsize=9)
    axs[0].legend(fontsize=7, ncol=2)
    w0s = sorted({r["w0"] for r in rows}, key=lambda s: -float(s.split("/")[1]))
    x = np.arange(len(w0s))
    wd = 0.8 / len(names)
    for i, k in enumerate(names):
        a = [min(r["alpha"] for r in rows if r["variant"] == k and r["w0"] == w) for w in w0s]
        axs[1].bar(x + (i - (len(names) - 1) / 2) * wd, a, width=wd * 0.9, color=COL[k], label=k, zorder=2)
    axs[1].set_xticks(x, w0s)
    axs[1].set_xlabel("imperfection amplitude w₀")
    axs[1].set_ylabel("knockdown α = λ_GNIA / λ_LBA (worse sign)")
    axs[1].legend(fontsize=8, ncol=len(names))
    fig.tight_layout()
    fig.savefig(FIG / "fig4_gnia.png")
    plt.close(fig)
    table("table_gnia.csv", ["variant", "w0", "sign", "angle_deg", "lambda_LBA", "lambda_GNIA", "alpha", "kind"],
          [[r["variant"], r["w0"], r["sign"], r["angle_deg"], round(r["lambda_LBA"], 2),
            round(r["lambda_GNIA"], 2), round(r["alpha"], 4), r["kind"]] for r in rows])


def fig_modes(h: float, lay: dict | None) -> None:
    p = RES / f"modes_h{h:g}.npz"
    if not p.exists():
        print(f"skip: {p.name} not found (run s07_exports.py first)")
        return
    z = np.load(p)
    names = [k for k in ("A", "C", "D", "E") if f"{k}_xyz" in z]
    fig, axs = plt.subplots(1, len(names), figsize=(2.6 * len(names) + 0.8, 3.0))
    for ax, k in zip(np.atleast_1d(axs), names):
        xyz, quads, amp = z[f"{k}_xyz"], z[f"{k}_quads"], z[f"{k}_amp"]
        polys = xyz[quads][:, :, :2]
        pc = PolyCollection(polys, array=amp[quads].mean(axis=1), cmap="Blues", edgecolor="none")
        pc.set_clim(0, 1)
        ax.add_collection(pc)
        if lay and k in lay:
            for ln in lay[k]:
                ln = np.asarray(ln)
                ax.plot(ln[:, 0], ln[:, 1], color=INK, lw=0.9)
        ax.set_aspect("equal")
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 12)
        ax.set_xticks([0, 6, 12])
        ax.set_yticks([0, 6, 12])
        ax.grid(False)
        ax.set_title(f"{k}: λ = {float(z[f'{k}_factor']):.1f} ({int(z[f'{k}_angle'])}°)", fontsize=9)
    fig.colorbar(pc, ax=np.atleast_1d(axs).tolist(), label="mode-1 |u| / max", shrink=0.8)
    fig.savefig(FIG / "fig5_modes.png", bbox_inches="tight")
    plt.close(fig)


def fig_angles(var: dict) -> None:
    v = var["variants"]
    names = [k for k in COL if k in v]
    ang = sorted(int(a) for a in v["A"]["lambda_c_by_angle"])
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for k in names:
        lam = [v[k]["lambda_c_by_angle"][str(a)] for a in ang]
        ax.plot(ang, lam, color=COL[k], marker=MRK[k], ms=6)
        ax.annotate(k, (ang[-1], lam[-1]), xytext=(6, 0), textcoords="offset points", va="center",
                    fontweight="bold")
    ax.set_xticks(ang)
    ax.set_xlabel("sway direction [deg from +x]")
    ax.set_ylabel("λ_c at V = V60")
    fig.tight_layout()
    fig.savefig(FIG / "fig6_by_angle.png")
    plt.close(fig)
    table("table_variants_V60.csv", ["variant", "lambda_g", "lambda_c", "worst_angle", "rib_length_m", "n_ribs",
                                     "crossings", "rib_depth_mm"],
          [[k, round(v[k]["lambda_g"], 2), round(v[k]["lambda_c"], 2), v[k]["worst_angle"],
            round(v[k].get("rib_length_m", 0), 1), v[k].get("n_ribs", 0), v[k].get("crossings", 0),
            round(v[k].get("rib_depth_mm", 0), 1)] for k in names])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.4)
    h = ap.parse_args().h
    FIG.mkdir(exist_ok=True)
    lay = load("layouts.json")
    if (sw := load(f"sweep_h{h:g}.json")):
        fig_sweep(sw)
    fig_E(h)
    if (g := load(f"gnia_h{h:g}.json")):
        fig_gnia(g)
    fig_modes(h, lay)
    if (var := load(f"variants_h{h:g}_r1.json")):
        fig_angles(var)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
