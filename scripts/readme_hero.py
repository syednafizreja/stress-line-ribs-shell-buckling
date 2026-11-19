"""Draw the README hero from the result files.

Left: linear buckling factor at the shared volume V60. Right: the governing
geometrically nonlinear factor. One PNG per GitHub theme, 1280 x 640 px, so
the light one is also the social preview. The plotted numbers are written to
docs/img/hero.json. Nothing here changes a result.

Colours are the variant identity colours in study/s06_figures.py, in study
order. They are not sorted by the result.

Run: python scripts/readme_hero.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "figures" / "table_variants_V60.csv"
E_JSON = ROOT / "results" / "E_best_d4_h0.4.json"
GNIA = ROOT / "results" / "gnia_h0.4.txt"
OUT = ROOT / "docs" / "img"

# Identity colours, fixed in study/s06_figures.py. Never re-ordered by rank.
COL = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a", "D": "#eda100",
       "Dd": "#e87ba4", "Dx": "#008300", "E": "#6250d6"}
LINEAR = ["A", "B", "C", "D", "Dd", "Dx", "E"]
NONLINEAR = ["A", "C", "D", "E"]
ROLE = {"A": "uniform", "B": "grid", "C": "gravity", "D": "envelope",
        "Dd": "diagonal", "Dx": "+x only", "E": "optimised"}

THEMES = {
    "light": {"bg": "#ffffff", "ink": "#0b0b0b", "muted": "#52514e", "faint": "#898781",
              "grid": "#e1e0d9", "spine": "#c3c2b7"},
    "dark": {"bg": "#0d1117", "ink": "#f0efec", "muted": "#c3c2b7", "faint": "#8f8d86",
             "grid": "#30363d", "spine": "#484f58"},
}


def load_linear() -> dict[str, float]:
    """lambda_c from the V60 table, plus E's objective from its result file."""
    with CSV.open(encoding="utf-8", newline="") as f:
        rows = {r["variant"]: float(r["lambda_c"]) for r in csv.DictReader(f)}
    e = json.loads(E_JSON.read_text(encoding="utf-8"))
    rows["E"] = float(e["objective"])
    missing = [k for k in LINEAR if k not in rows]
    if missing:
        raise SystemExit(f"missing variants: {missing}")
    return rows


def load_gnia() -> dict[str, dict[str, str]]:
    """Minimum over amplitude of the weaker-sign row in the GNIA record."""
    lines = GNIA.read_text(encoding="utf-8").splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("governing"))
    best: dict[str, dict[str, str]] = {}
    for ln in lines[start + 1:]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        var, w0, lam, alpha = parts[:4]
        if var not in best or float(lam) < float(best[var]["lam"]):
            best[var] = {"w0": w0, "lam": lam, "alpha": alpha}
    missing = [k for k in NONLINEAR if k not in best]
    if missing:
        raise SystemExit(f"GNIA record has no governing row for {missing}")
    return best


def derived(lam: dict[str, float], gnia: dict[str, dict[str, str]]) -> dict:
    """Comparisons a reader can recompute. Stored so the page and the figure match."""
    a, e = lam["A"], lam["E"]
    gain = e - a

    def pct(v: float) -> float:
        return 100.0 * (v - a) / a

    def recover(v: float) -> float:
        return 100.0 * (v - a) / gain

    g = {k: float(gnia[k]["lam"]) for k in NONLINEAR}
    return {
        "lambda_c": {k: lam[k] for k in LINEAR},
        "gain_E_minus_A": gain,
        "percent_vs_A": {k: pct(lam[k]) for k in LINEAR if k != "A"},
        "percent_of_E_gain": {k: recover(lam[k]) for k in ("B", "C", "D", "Dd", "Dx")},
        "gnia": {k: {"w0": gnia[k]["w0"], "lambda": g[k], "alpha": gnia[k]["alpha"]} for k in NONLINEAR},
        "gnia_E_vs_A_percent": 100.0 * (g["E"] - g["A"]) / g["A"],
    }


def _fit(fig, text: str, x: float, y: float, color: str, weight: str, start: float, floor: float) -> None:
    t = fig.text(x, y, text, fontsize=start, fontweight=weight, color=color, ha="left", va="top")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    size = start
    while size > floor:
        bb = t.get_window_extent(renderer)
        if bb.x1 <= fig.bbox.xmax - 28:
            return
        size -= 0.5
        t.set_fontsize(size)
        fig.canvas.draw()


def fonts() -> list[str]:
    """Arial, then its metric twin, then DejaVu. Skip names that are not installed."""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    names = [n for n in ("Arial", "Liberation Sans", "DejaVu Sans") if n in have]
    return names or ["DejaVu Sans"]


def draw(lam: dict[str, float], gnia: dict[str, dict[str, str]], theme: str) -> Path:
    th = THEMES[theme]
    plt.rcParams.update({"font.family": fonts(), "font.size": 11, "axes.unicode_minus": False})
    fig = plt.figure(figsize=(12.8, 6.4), dpi=100, facecolor=th["bg"])
    headline = (f"Stress-line ribs fall below the plain shell: {lam['D']:.2f} against {lam['A']:.2f}. "
                f"Optimised thickness reaches {lam['E']:.4f}.")
    _fit(fig, headline, 0.045, 0.955, th["ink"], "bold", 16.5, 13)
    fig.text(0.045, 0.875,
             "Same concrete volume as a uniform 60 mm shell.  Left: linear factor, worst direction.  "
             "Right: governing nonlinear factor, weaker imperfection sign, same scale.",
             fontsize=11, color=th["muted"], ha="left", va="top")

    ymax = max(lam.values()) * 1.22
    panels = [
        (0.055, LINEAR, {k: lam[k] for k in LINEAR}, "Linear buckling factor", False),
        (0.545, NONLINEAR, {k: float(gnia[k]["lam"]) for k in NONLINEAR},
         "Nonlinear, governing value, same scale", True),
    ]
    for x0, keys, values, title, nonlinear in panels:
        ax = fig.add_axes([x0, 0.175, 0.40, 0.56], facecolor=th["bg"])
        xpos = np.arange(len(keys))
        heights = [values[k] for k in keys]
        ax.bar(xpos, heights, width=0.68, color=[COL[k] for k in keys], zorder=3, linewidth=0)
        ax.axhline(values["A"], color=th["muted"], lw=0.8, ls=(0, (3, 2)), zorder=4)
        for x, k in zip(xpos, keys):
            label = f"{values[k]:.4f}" if k == "E" and not nonlinear else f"{values[k]:.2f}"
            ax.text(x, values[k] + ymax * 0.02, label, ha="center", va="bottom",
                    fontsize=9.5, color=th["ink"], zorder=5)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.55, len(keys) - 0.45)
        ax.set_xticks(xpos, [f"{k}\n{ROLE[k]}" for k in keys])
        ax.tick_params(axis="x", length=0, pad=4, colors=th["muted"], labelsize=9)
        ax.tick_params(axis="y", length=3, direction="out", colors=th["muted"], labelsize=9)
        ax.yaxis.grid(True, color=th["grid"], lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for name, spine in ax.spines.items():
            spine.set_visible(name in ("left", "bottom"))
            spine.set_color(th["spine"])
            spine.set_linewidth(0.6)
        ax.set_title(title, loc="left", color=th["ink"], fontsize=12.5, fontweight="bold", pad=8)
        if x0 < 0.1:
            ax.set_ylabel("Factor on self-weight + 0.2 g", color=th["muted"], fontsize=10)
        else:
            ax.set_yticklabels([])
    amps = ", ".join(f"{k} {gnia[k]['w0']}" for k in NONLINEAR)
    lines = [
        "Bars start at zero and stay in study order. Dashed line: design A on that panel. The number is the file value.",
        f"Nonlinear amplitudes: {amps}. Data: figures/table_variants_V60.csv, results/E_best_d4_h0.4.json, results/gnia_h0.4.txt.",
    ]
    for y, line in ((0.058, lines[0]), (0.026, lines[1])):
        _fit(fig, line, 0.045, y, th["faint"], "normal", 8.5, 7.5)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"hero_{theme}.png"
    fig.savefig(path, dpi=100, facecolor=th["bg"], metadata={"Software": None})
    plt.close(fig)
    return path


def main() -> None:
    lam = load_linear()
    gnia = load_gnia()
    for theme in THEMES:
        print(draw(lam, gnia, theme))
    payload = {
        "sources": [str(CSV.relative_to(ROOT)).replace("\\", "/"),
                    str(E_JSON.relative_to(ROOT)).replace("\\", "/"),
                    str(GNIA.relative_to(ROOT)).replace("\\", "/")],
        "rule_linear": "lambda_c from the V60 table; E is the objective in E_best_d4_h0.4.json",
        "rule_nonlinear": "minimum over imperfection amplitude of the weaker-sign row",
        "order_linear": LINEAR,
        "order_nonlinear": NONLINEAR,
        **derived(lam, gnia),
    }
    (OUT / "hero.json").write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    from PIL import Image
    for theme in THEMES:
        im = Image.open(OUT / f"hero_{theme}.png")
        if im.size != (1280, 640):
            raise SystemExit(f"hero_{theme}.png is {im.size}, expected 1280x640")
    print("wrote", OUT / "hero.json")


if __name__ == "__main__":
    main()
