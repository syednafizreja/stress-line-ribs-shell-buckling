"""Do the ribs work, or is the ribbed shell governed by its panels?

At equal volume the ribbed variants buckle at a lower lambda_c than the uniform shell A, and D
lands close to a bare 50 mm shell carrying the V60 weight. That can be real (the thin panels
between ribs buckle locally, with the ribs as near-nodal lines) or a fault (the B32R ribs are
not coupled to the S8R shell and add mass but no stiffness). Three discriminating checks, h = 0.6:

  1. bare 50 mm shell with density scaled by 60/50, so it carries exactly the V60 weight;
  2. mode-1 shape of C and D: mean and peak |u| on the rib nodes relative to the whole shell
     (ribs that work are close to nodal lines of a panel mode);
  3. C and D with three times the rib depth (volume not held): working ribs raise lambda_c only
     if the governing mode involves them; a pure panel mode barely changes.

Usage: python verify/rib_effectiveness.py   (writes verify/out/rib_effectiveness.txt and mode npz)
"""
from __future__ import annotations

import dataclasses
import json
import tempfile
from pathlib import Path

import numpy as np

from ribshell.constants import RHO, T_BASE_RIBBED, T_REF
from ribshell.fe import Model, Rib
from ribshell.form import shell_surface
from ribshell.gnia import lba_mode
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).parent / "out"
H = 0.6


def rib_share(m: Model, mode: dict[int, np.ndarray]) -> tuple[float, float]:
    rib_nodes = {n for r in m.ribs for tri in m.mesh.ribs[r.rib_id] for n in tri}
    free = [n for n in m.mesh.nodes if n not in m.mesh.edge_nodes]
    amp = {n: float(np.linalg.norm(mode.get(n, np.zeros(3)))) for n in free}
    peak = max(amp.values())
    on = [amp[n] for n in rib_nodes if n in amp]
    return float(np.mean(on) / np.mean(list(amp.values()))), float(max(on) / peak)


def main() -> None:
    s = shell_surface()
    mesh = build_mesh(s, H)
    v60 = mesh.area() * T_REF
    lay = json.load(open(ROOT / "results" / "layouts.json"))
    rows = [f"rib effectiveness, h = {H} m, G + 0.2 g at 0 deg, lowest positive factor", ""]
    modes = {}
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        a = lba_mode(Model(mesh, T_REF, surface=s), wd, "A", 0)
        bare = lba_mode(Model(mesh, T_BASE_RIBBED, surface=s, rho=RHO * T_REF / T_BASE_RIBBED), wd, "bare", 0)
        rows += [f"A uniform 60 mm                          lambda_c0 = {a.factor:8.2f}",
                 f"bare 50 mm shell carrying the V60 weight lambda_c0 = {bare.factor:8.2f}", ""]
        modes["A"] = a
        for n in ("C", "D"):
            m, met = ribbed_model(Layout(n, [np.array(x) for x in lay[n]]), s, v60, h=H)
            r1 = lba_mode(m, wd, n, 0)
            mean_s, peak_s = rib_share(m, r1.mode)
            deep = dataclasses.replace(m, ribs=[Rib(r.rib_id, r.width, 3 * r.depth) for r in m.ribs])
            r3 = lba_mode(deep, wd, n + "3", 0)
            modes[n] = r1
            rows += [f"{n}: {met['n_ribs']} ribs, {met['rib_length_m']:.1f} m, depth {met['rib_depth_mm']:.0f} mm",
                     f"   lambda_c0 = {r1.factor:8.2f}   (x {r1.factor / bare.factor:.3f} of the bare shell)",
                     f"   mode 1 on rib nodes: mean |u| / shell mean = {mean_s:.3f}, max |u| / shell peak = {peak_s:.3f}",
                     f"   ribs 3x deeper ({3 * met['rib_depth_mm']:.0f} mm, heavier): lambda_c0 = {r3.factor:8.2f}"
                     f"  (x {r3.factor / r1.factor:.3f})"]
    OUT.mkdir(exist_ok=True)
    (OUT / "rib_effectiveness.txt").write_text("\n".join(rows) + "\n")
    np.savez_compressed(OUT / "rib_effectiveness_modes.npz",
                        **{f"{k}_ids": np.array(list(v.mode)) for k, v in modes.items()},
                        **{f"{k}_u": np.array(list(v.mode.values())) for k, v in modes.items()})
    print("\n".join(rows))


if __name__ == "__main__":
    main()
