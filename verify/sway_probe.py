"""Does lateral sway matter for buckling?  Edge-pinned vs corner-supported funicular shell.

For each support case and each a_h (0.15, 0.2, 0.3, 0.4 g) the script reports lambda_g (gravity
alone), lambda_c (G + a_h together, worst of 0 and 45 degrees: the uniform square is D4-symmetric,
so these two directions cover the envelope) and the sway ratio lambda_c / lambda_g.

Corner-supported case (brief section 10): fixed-plan funicular form with edge arches rising 0.5 x
the crown, a 250 x 400 mm edge beam hanging below the free edges (it carries the edge thrust that a
fixed-plan form cannot balance in plan), and pinned 0.30 m x 0.30 m column pads at the corners.

Decision rule (STATUS, 17 Nov): lambda_c/lambda_g <= ~0.7 at 0.2 g makes the corner case a core case.

Usage:  python verify/sway_probe.py      (writes verify/out/sway_probe.txt)
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from ribshell.constants import SPAN, T_REF
from ribshell.fe import Model, Rib, run_buckle
from ribshell.form import corner_surface, shell_surface
from ribshell.mesh import build_mesh

OUT = Path(__file__).parent / "out"
H = 0.4
A_HS = (0.15, 0.20, 0.30, 0.40)
EDGE_BEAM = (0.25, 0.40)   # width, depth [m]
PAD = 0.30                 # corner column pad side [m]


def edge_pinned() -> Model:
    s = shell_surface()
    return Model(build_mesh(s, H), T_REF, surface=s)


def corner_supported() -> Model:
    s = corner_surface()
    L = SPAN
    edges = [[(0, 0), (L, 0)], [(L, 0), (L, L)], [(L, L), (0, L)], [(0, L), (0, 0)]]
    import numpy as np
    mesh = build_mesh(s, H, ribs=[np.array(e, float) for e in edges])
    sup = {n for n, (x, y, _) in mesh.nodes.items()
           if min(x, L - x) <= PAD + 1e-9 and min(y, L - y) <= PAD + 1e-9
           and (min(x, L - x) < 1e-6 * L or min(y, L - y) < 1e-6 * L)}
    ribs = [Rib(i, *EDGE_BEAM) for i in range(4)]
    return Model(mesh, T_REF, ribs=ribs, surface=s, supports=sup)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    lines = [f"Sway probe, t = {T_REF*1e3:.0f} mm, h = {H} m, lambda_1 of 10 modes",
             "case            a_h   lam_g     lam_c(0)  lam_c(45)  lam_c    lam_c/lam_g"]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        for name, build in (("edge-pinned", edge_pinned), ("corner+edgebeam", corner_supported)):
            m = build()
            lg = run_buckle(m, d, f"{name[:4]}_g", "gravity").factors[0]
            for a in A_HS:
                lc = [run_buckle(m, d, f"{name[:4]}_{int(a*100)}_{ang}", "combined", ang, a_h=a).factors[0]
                      for ang in (0, 45)]
                lines.append(f"{name:15s} {a:4.2f}  {lg:8.2f}  {lc[0]:8.2f}  {lc[1]:8.2f}  "
                             f"{min(lc):8.2f}  {min(lc)/lg:6.3f}")
                print(lines[-1], flush=True)
    (OUT / "sway_probe.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
