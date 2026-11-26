"""Publication figures and booktabs tables from paper/FACTS.md.

Numbers are the verified sheet (commit 2025f66). Equilibrium histories,
the mesh series, and the D4 thickness grid are read from the files that
sheet names. Missing nonlinear runs stay missing.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "paper" / "figures"
PRE = FIG / "preview"
TAB = ROOT / "paper" / "tables"

NAVY = "#1f3a5f"
SLATE = "#4c5d73"
BRICK = "#a34b3c"
TEAL = "#2a6b66"
GOLD = "#b0893f"
GRAY = "#6b7280"
GRID = "#e6e8eb"
INK = "#1c1c1c"

# Eight-direction linear factors, section 3. Baseline A = 341.2093.
LINEAR = {
    "A": 341.2093,
    "B": 395.9731,
    "C": 249.1289,
    "D": 235.0948,
    "Dd": 268.9580,
    "Dx": 193.1005,
    "E": 429.0206,
}
# Weaker-sign minimum over L/300, L/500, L/1000. B, Dd, Dx were not run.
GNIA_MIN = {"A": 118.202975, "C": 91.295285, "D": 89.771960, "E": 135.958464}
# Peak-down alphas at the three amplitudes (all limit points).
ALPHA = {
    "A": {"L/1000": 0.434151, "L/500": 0.346424, "L/300": 0.356874},
    "C": {"L/1000": 0.502097, "L/500": 0.366453, "L/300": 0.423750},
    "D": {"L/1000": 0.503467, "L/500": 0.418063, "L/300": 0.381826},
    "E": {"L/1000": 0.386781, "L/500": 0.316904, "L/300": 0.327269},
}
ORDER = ["A", "B", "C", "D", "Dd", "Dx", "E"]
NAMES = {
    "A": "A  plain 60 mm",
    "B": "B  orthogonal grid",
    "C": "C  gravity stress lines",
    "D": "D  moment stress lines",
    "Dd": "Dd diagonal moment",
    "Dx": "Dx +x moment lines",
    "E": "E  optimised thickness",
}
PARAMS = {
    "A": "uniform 60 mm",
    "B": "14 ribs, 62.4 mm, 49 crossings",
    "C": "8 ribs, 157.6 mm, 0 crossings",
    "D": "12 ribs, 208.8 mm, 0 crossings",
    "Dd": "8 ribs, 400.2 mm, 0 crossings",
    "Dx": "5 ribs, 505.1 mm, 0 crossings",
    "E": "40.0--80.6 mm, 5 x 5 zones",
}

# h = 0.6 m, theta = 0 only. Section 4.
RIB_CHECK = [
    ("A, 60 mm", 341.41),
    ("Bare 50 mm", 221.23),
    ("C, 158 mm", 249.75),
    ("D, 209 mm", 237.68),
    ("C, depth x3", 227.56),
    ("D, depth x3", 179.92),
]

SCORDELIS_N = np.array([8, 16, 24, 32])
SCORDELIS_U = np.array([0.30325, 0.30196, 0.30192, 0.30192])
SCORDELIS_REF = 0.3024
ZOELLY_REF = 3.8891
ZOELLY_FE = 3.9224

MESH_H = np.array([0.6, 0.4, 0.3, 0.2])
MESH = {
    "A": np.array([341.41, 341.21, 341.18, 341.16]),
    "C": np.array([249.75, 249.13, 249.05, 248.99]),
    "D": np.array([237.68, 235.11, 234.89, 234.69]),
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.edgecolor": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )


def tidy(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="both", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    PRE.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(PRE / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def pct(value: float, base: float) -> float:
    return 100.0 * (value - base) / base


def fig1() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))

    ax = axes[0]
    ax.axhline(SCORDELIS_REF, color=GRAY, lw=0.8, ls="--", zorder=1)
    ax.plot(SCORDELIS_N, SCORDELIS_U, color=NAVY, marker="o", ms=4.5, lw=1.1, zorder=2)
    ax.set_xlabel("Elements along each edge")
    ax.set_ylabel("Free-edge vertical displacement")
    ax.set_title("Scordelis–Lo roof, S8R")
    ax.set_xticks(SCORDELIS_N)
    ax.set_ylim(0.30155, 0.30355)
    err = 100.0 * (SCORDELIS_U[-1] - SCORDELIS_REF) / SCORDELIS_REF
    ax.text(
        0.98,
        0.96,
        f"32×32: {SCORDELIS_U[-1]:.5f}\n{err:.2f}% vs 0.3024",
        transform=ax.transAxes,
        fontsize=7,
        color=SLATE,
        ha="right",
        va="top",
    )
    tidy(ax)

    ax = axes[1]
    labels = ["Zoelly\nformula", "S8R\nhemisphere"]
    vals = [ZOELLY_REF, ZOELLY_FE]
    bars = ax.bar(labels, vals, color=[GRAY, TEAL], width=0.55, zorder=2)
    ax.set_ylabel("Critical pressure (MPa)")
    ax.set_title("Zoelly spherical cap")
    ax.set_ylim(0, 6.4)
    for bar, val in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + 0.08,
            f"{val:.4f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    err_z = 100.0 * (ZOELLY_FE - ZOELLY_REF) / ZOELLY_REF
    ax.text(
        0.5,
        0.90,
        f"R = 5 m, t = 50 mm, 27 972 S8R, {err_z:+.2f}%",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=6.5,
        color=SLATE,
    )
    tidy(ax)

    ax = axes[2]
    series = (("A", NAVY), ("C", TEAL), ("D", BRICK))
    for key, color in series:
        ax.plot(MESH_H, MESH[key], color=color, marker="o", ms=3.5, lw=1.1, label=key, zorder=2)
    ax.set_xlabel("Mesh size h (m)")
    ax.set_ylabel(r"$\lambda$ at $0^\circ$")
    ax.set_title("Study mesh")
    ax.invert_xaxis()
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)
    fig.tight_layout(w_pad=1.4)
    save(fig, "fig1_benchmarks")


def fig2() -> None:
    fig, ax = plt.subplots(figsize=(7.1, 3.55))
    ids = list(reversed(ORDER))
    y = np.arange(len(ids))
    h = 0.36
    base = LINEAR["A"]
    gbase = GNIA_MIN["A"]
    lin = [LINEAR[k] for k in ids]
    ax.barh(y + h / 2, lin, height=h, color=NAVY, zorder=2, label=r"Linear $\lambda_c$")
    gnia_y = []
    gnia_v = []
    for i, k in enumerate(ids):
        if k in GNIA_MIN:
            gnia_y.append(y[i] - h / 2)
            gnia_v.append(GNIA_MIN[k])
        else:
            ax.text(6, y[i] - h / 2, "not run", va="center", fontsize=7, color=GRAY)
    ax.barh(
        gnia_y,
        gnia_v,
        height=h,
        color=BRICK,
        zorder=2,
        label=r"Mode-1 $\lambda_{\mathrm{GNIA}}$ (weaker sign, min. amplitude)",
    )
    for i, k in enumerate(ids):
        shift = pct(LINEAR[k], base)
        tag = "baseline" if k == "A" else f"{shift:+.2f}%"
        ax.text(LINEAR[k] + 6, y[i] + h / 2, tag, va="center", fontsize=7, color=NAVY)
        if k in GNIA_MIN:
            gshift = pct(GNIA_MIN[k], gbase)
            gtag = "baseline" if k == "A" else f"{gshift:+.2f}%"
            ax.text(
                GNIA_MIN[k] + 6,
                y[i] - h / 2,
                gtag,
                va="center",
                fontsize=7,
                color=BRICK,
            )
    ax.set_yticks(y)
    ax.set_yticklabels([NAMES[k] for k in ids])
    ax.set_xlabel("Load factor")
    ax.set_xlim(0, 540)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.42, -0.18), ncol=1)
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig2_layout_buckling")


def _histories():
    path = ROOT / "results" / "gnia_h0.4.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    return data["rows"]


def fig3() -> None:
    rows = _histories()
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.15))
    colors = {"A": NAVY, "C": TEAL, "D": BRICK, "E": GOLD}
    labels = {
        "A": "A plain shell",
        "C": "C stress-line ribs",
        "D": "D stress-line ribs",
        "E": "E optimised thickness",
    }
    ax = axes[0]
    for row in rows:
        if row["w0"] != "L/500":
            continue
        if row["variant"] not in colors:
            continue
        hist = np.asarray(row["history"], dtype=float)
        lam = hist[:, 0]
        umm = hist[:, 1] * 1000.0
        solid = row["kind"] == "limit point"
        ax.plot(
            umm,
            lam,
            color=colors[row["variant"]],
            lw=1.15,
            ls="-" if solid else "--",
            zorder=2,
        )
    handles = [
        plt.Line2D([0], [0], color=colors[k], lw=1.15, label=labels[k])
        for k in ("A", "C", "D", "E")
    ]
    handles.append(plt.Line2D([0], [0], color=INK, lw=1.0, label="limit point"))
    handles.append(
        plt.Line2D([0], [0], color=INK, lw=1.0, ls="--", label="solver stop")
    )
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=7)
    ax.set_xlabel(r"Displacement at the mode-1 peak node (mm)")
    ax.set_ylabel(r"Load factor $\lambda$")
    ax.set_title(r"Equilibrium paths at $w_0 = L/500$")
    tidy(ax)

    ax = axes[1]
    amps = ["L/1000", "L/500", "L/300"]
    x = np.array([1 / 1000, 1 / 500, 1 / 300])
    for k in ("A", "C", "D", "E"):
        y = [ALPHA[k][a] for a in amps]
        ax.plot(x, y, color=colors[k], marker="o", ms=4, lw=1.15, label=labels[k], zorder=2)
    ax.set_xlabel(r"Imperfection amplitude $w_0 / L$")
    ax.set_ylabel(r"Knockdown $\alpha = \lambda_{\mathrm{GNIA}} / \lambda_{\mathrm{LBA}}$")
    ax.set_title("Weaker sign (peak-down)")
    ax.set_xticks(x)
    ax.set_xticklabels(["1/1000", "1/500", "1/300"])
    ax.set_ylim(0.25, 0.58)
    ax.legend(frameon=False, fontsize=7)
    tidy(ax)
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig3_imperfection_sensitivity")


def fig4() -> None:
    raw = json.loads(
        (ROOT / "results" / "E_best_d4_h0.4.json").read_text(encoding="utf-8")
    )
    zones = np.asarray(raw["zone_thickness_m"], dtype=float).reshape(5, 5) * 1000.0

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.25), gridspec_kw={"width_ratios": [1.15, 1]})

    ax = axes[0]
    labels = [a[0] for a in RIB_CHECK]
    vals = [a[1] for a in RIB_CHECK]
    colors = [NAVY, GRAY, TEAL, BRICK, TEAL, BRICK]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color=colors, height=0.62, zorder=2)
    for yi, val in zip(y, vals):
        ax.text(val + 4, yi, f"{val:.2f}", va="center", fontsize=7, color=SLATE)
    ax.axvline(221.23, color=GRAY, lw=0.6, ls=":", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$\lambda$ at $0^\circ$, $h = 0.6$ m")
    ax.set_xlim(0, 430)
    ax.set_title("Ribs on a 50 mm shell")
    tidy(ax)

    ax = axes[1]
    # ix * 5 + iy; rows of the image follow y.
    image = ax.imshow(
        zones.T,
        origin="lower",
        cmap="cividis",
        vmin=40,
        vmax=85,
        extent=(0, 12, 0, 12),
        zorder=2,
    )
    for ix in range(5):
        for iy in range(5):
            ax.text(
                1.2 + 2.4 * ix,
                1.2 + 2.4 * iy,
                f"{zones[ix, iy]:.0f}",
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if zones[ix, iy] < 70 else INK,
            )
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.ax.tick_params(labelsize=7, width=0.6, length=2)
    cbar.set_label("Thickness (mm)", fontsize=8)
    cbar.outline.set_linewidth(0.6)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("E, same volume")
    ax.set_xticks([0, 6, 12])
    ax.set_yticks([0, 6, 12])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.tight_layout(w_pad=1.6)
    save(fig, "fig4_thickness_stiffness_mechanism")


def write_tables() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    body = [
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Check & Mesh & Reference & Finite element & Error \\",
        r"\midrule",
        r"Scordelis--Lo, free-edge $u$ & $8\times 8$ & 0.3024 & 0.30325 & $+0.28\%$ \\",
        r" & $16\times 16$ & 0.3024 & 0.30196 & $-0.15\%$ \\",
        r" & $24\times 24$ & 0.3024 & 0.30192 & $-0.16\%$ \\",
        r" & $32\times 32$ & 0.3024 & 0.30192 & $-0.16\%$ \\",
        r"Zoelly cap, $p_{\mathrm{cr}}$ (MPa) & 27\,972 S8R & 3.8891 & 3.9224 & $+0.86\%$ \\",
        r"\midrule",
        r"\multicolumn{5}{l}{Study mesh, $\lambda$ at $0^\circ$, lowest positive of 10 modes} \\",
        r"\midrule",
        r"$h$ (m) & A $\lambda_{c0}$ & C $\lambda_{c0}$ & D $\lambda_{c0}$ & \\",
    ]
    for h, a, c, d in zip(MESH_H, MESH["A"], MESH["C"], MESH["D"]):
        body.append(f"{h:.1f} & {a:.2f} & {c:.2f} & {d:.2f} & \\\\")
    body += [
        r"\bottomrule",
        r"\end{tabular}",
        "",
    ]
    (TAB / "tab1_benchmarks.tex").write_text("\n".join(body), encoding="utf-8")

    base = LINEAR["A"]
    gbase = GNIA_MIN["A"]
    lines = [
        r"\begin{tabular}{lrrrrl}",
        r"\toprule",
        r"Layout & $\lambda_c$ & vs.\ A & $\lambda_{\mathrm{GNIA}}$ & vs.\ A & Governing geometry \\",
        r"\midrule",
    ]
    for k in ORDER:
        lin = LINEAR[k]
        lp = pct(lin, base)
        ltag = "baseline" if k == "A" else f"{lp:+.2f}\\%"
        if k in GNIA_MIN:
            g = GNIA_MIN[k]
            gp = pct(g, gbase)
            gtag = "baseline" if k == "A" else f"{gp:+.2f}\\%"
            gcell = f"{g:.2f}"
        else:
            gcell = "---"
            gtag = "not run"
        lines.append(
            f"{k} & {lin:.2f} & {ltag} & {gcell} & {gtag} & {PARAMS[k]} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    (TAB / "tab2_layouts.tex").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    style()
    fig1()
    fig2()
    fig3()
    fig4()
    write_tables()
    names = [
        "paper/figures/fig1_benchmarks.pdf",
        "paper/figures/fig2_layout_buckling.pdf",
        "paper/figures/fig3_imperfection_sensitivity.pdf",
        "paper/figures/fig4_thickness_stiffness_mechanism.pdf",
        "paper/figures/preview/fig1_benchmarks.png",
        "paper/figures/preview/fig2_layout_buckling.png",
        "paper/figures/preview/fig3_imperfection_sensitivity.png",
        "paper/figures/preview/fig4_thickness_stiffness_mechanism.png",
        "paper/tables/tab1_benchmarks.tex",
        "paper/tables/tab2_layouts.tex",
    ]
    print("generated:")
    for name in names:
        path = ROOT / name
        print(f"  {name}  {path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
