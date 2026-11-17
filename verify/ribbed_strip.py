"""Check the rib model (S8R shell + B32R beam with OFFSET2 below the mid-surface) on a ribbed slab.

Flat 6 m x 6 m slab, t = 50 mm, eight ribs 150 x 250 mm along x at y = 0.375 + 0.75 k, so each rib
with its 0.75 m flange strip is one T-beam. Roller supports (u_z = 0) on x = 0 and x = a, free
edges at y = 0, a, self-weight. A pinned in-plane support would add a spurious end moment from the
eccentric rib, so u_x, u_y are restrained at two points only.

Reference: simply supported T-beam, flange b x t on the mid-surface, web w x d hanging below,
    w_mid = 5 q L^4 / (384 E I) + q L^2 / (8 kappa G A_web),   kappa = 5/6,
with q the self-weight per strip. Plate action of the flange (1 - nu^2) and residual shear lag
make the FE slightly stiffer or softer; pass band 5 % (scope section 7).

Usage: python verify/ribbed_strip.py   (writes verify/out/ribbed_strip.txt)
"""
from __future__ import annotations

import re
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

from ribshell.constants import E, G, NU, RHO
from ribshell.fe import Model, Rib, _grav, _model_cards, _set
from ribshell.form import Surface
from ribshell.mesh import build_mesh

OUT = Path(__file__).parent / "out"
A, T, W, D, NR, H = 6.0, 0.05, 0.15, 0.25, 8, 0.1


def flat(a: float) -> Surface:
    g = np.linspace(0, a, 5)
    return Surface(None, RectBivariateSpline(g, g, np.zeros((5, 5)), kx=3, ky=3), a)


def reference() -> tuple[float, dict]:
    b = A / NR
    af, aw = b * T, W * D
    zf, zw = 0.0, -(T / 2 + D / 2)
    zc = (af * zf + aw * zw) / (af + aw)
    i = b * T**3 / 12 + af * (zf - zc) ** 2 + W * D**3 / 12 + aw * (zw - zc) ** 2
    q = RHO * G * (af + aw)
    g_mod = E / (2 * (1 + NU))
    wb = 5 * q * A**4 / (384 * E * i)
    ws = q * A**2 / (8 * (5 / 6) * g_mod * aw)
    return wb + ws, {"I_m4": i, "q_N_per_m": q, "w_bend_mm": wb * 1e3, "w_shear_mm": ws * 1e3, "z_na_mm": zc * 1e3}


def _solve(mesh, m, cards, ys, tag):
    from cecore.env import run_ccx

    k = cards.index("EDGE,1,3")
    cards[k] = "EDGE,3,3"
    near = lambda p: min(mesh.nodes, key=lambda n: np.hypot(*(np.array(mesh.nodes[n][:2]) - p)))  # noqa: E731
    n1, n2 = near((0.0, A / 2)), near((A, A / 2))
    cards += [f"{n1},1,2", f"{n2},2,2"]
    mids = [near((A / 2, y)) for y in ys]
    deck = cards + _set("MID", mids) + ["*STEP", "*STATIC"] + \
        [ln for ln in _grav(G, (0.0, 0.0, -1.0)) if m.ribs or "RIBS" not in ln] + \
        ["*NODE PRINT,NSET=MID", "U", "*END STEP"]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / f"{tag}.inp").write_text("\n".join(deck) + "\n")
        r = run_ccx(tag, d, threads=4)
        assert r.returncode == 0, (r.stdout + r.stderr)[-800:]
        txt = (d / f"{tag}.dat").read_text()
    uz = {int(a): float(c) for a, c in re.findall(r"^\s*(\d+)\s+\S+\s+\S+\s+(\S+)\s*$", txt, re.M)}
    return np.array([-uz[n] for n in mids])


def beam_ribs(s, ys):
    """S8R slab + B32R ribs with OFFSET2 (the model used in the study)."""
    lines = [np.array([[0.0, y], [A, y]]) for y in ys]
    mesh = build_mesh(s, H, ribs=lines)
    sup = {n for n, (x, y, _) in mesh.nodes.items() if x < 1e-9 or x > A - 1e-9}
    m = Model(mesh, T, ribs=[Rib(i, W, D) for i in range(NR)], surface=s, supports=sup)
    return _solve(mesh, m, _model_cards(m), ys, "beam"), f"{len(mesh.quads)} S8R + " \
        f"{sum(len(mesh.ribs[i]) for i in range(NR))} B32R"


def shell_strips(s, ys):
    """Ribs as thickened S8R strips (t + d over the rib width), top flush with the slab.

    The strip edges y_r +- w/2 are embedded so the strip is meshed exactly. The strip section is
    offset so its mid-surface lies (t + d)/2 - t/2 below the nodes: OFFSET = -(d/2)/(t + d) times
    the sign of the element normal (CalculiX/Abaqus convention: offset of the mid-surface from the
    reference surface, in thicknesses, along the normal). The flipped T has the same
    I, so a sign error could not show in this check; the sign is set from the normal anyway.
    """
    lines = [np.array([[0.0, y + sgn * W / 2], [A, y + sgn * W / 2]]) for y in ys for sgn in (-1, 1)]
    hs = W / 2
    mesh = build_mesh(s, hs, ribs=lines)
    sup = {n for n, (x, y, _) in mesh.nodes.items() if x < 1e-9 or x > A - 1e-9}
    cy = {e: np.mean([mesh.nodes[n][1] for n in c[:4]]) for e, c in mesh.quads.items()}
    strip = {e for e, y in cy.items() if min(abs(y - ys)) < W / 2}
    thick = {e: (T + D if e in strip else T) for e in mesh.quads}
    m = Model(mesh, thick, surface=s, supports=sup)
    cards = _model_cards(m)
    e0 = next(iter(strip))
    p = [np.array(mesh.nodes[n]) for n in mesh.quads[e0][:3]]
    nz = np.sign(np.cross(p[1] - p[0], p[2] - p[0])[2])
    off = -nz * (D / 2) / (T + D)
    for i, ln in enumerate(cards):
        if ln.startswith("*SHELL SECTION") and abs(float(cards[i + 1]) - (T + D)) < 1e-6:
            cards[i] = ln + f",OFFSET={off:.6f}"
    return _solve(mesh, m, cards, ys, "strips"), f"{len(mesh.quads)} S8R ({len(strip)} in strips), h={hs} m"


def main() -> None:
    s = flat(A)
    ys = 0.5 * A / NR + np.arange(NR) * A / NR
    w_ref, info = reference()
    out = [f"ribbed slab a={A} m t={T*1e3:.0f} mm, {NR} ribs {W*1e3:.0f}x{D*1e3:.0f} mm, roller supports",
           "reference T-beam: " + ", ".join(f"{k}={v:.4g}" for k, v in info.items()),
           f"w_ref = {w_ref*1e3:.4f} mm"]
    ok = True
    inner = {}
    for name, fn in (("B32R beam ribs (study model)", beam_ribs), ("thickened S8R strips", shell_strips)):
        w_fe, desc = fn(s, ys)
        inner[name] = w_fe[2:-2].mean()
        ratio = inner[name] / w_ref
        ok &= abs(ratio - 1) < 0.05
        out += [f"{name}: {desc}",
                "  FE midspan w per rib [mm]: " + " ".join(f"{v*1e3:.4f}" for v in w_fe),
                f"  inner four ribs = {inner[name]*1e3:.4f} mm, ratio FE/ref = {ratio:.4f}"]
    a, b = inner.values()
    out.append(f"beam / strip = {a / b:.4f}")
    out.append("PASS" if ok else "FAIL")
    print("\n".join(out))
    OUT.mkdir(exist_ok=True)
    (OUT / "ribbed_strip.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
