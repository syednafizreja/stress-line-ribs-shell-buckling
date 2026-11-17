"""Beam-shell coupling on the curved shell: a point load on a rib node, with ribs of growing depth.

Layout D at h = 0.6, no self-weight, 10 kN down at the rib node nearest the plan centre of the
longest rib. If the B32R ribs are coupled to the S8R shell (offset, orientation), the deflection
under the load falls steeply as the rib deepens (EI grows with d^3); an uncoupled rib would leave
it at the bare-shell value.

Usage: python verify/rib_coupling_curved.py   (writes verify/out/rib_coupling_curved.txt)
"""
from __future__ import annotations

import dataclasses
import json
import re
import tempfile
from pathlib import Path

import numpy as np

from cecore.env import run_ccx
from ribshell.constants import T_REF
from ribshell.fe import Rib, _model_cards
from ribshell.form import shell_surface
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).parent / "out"


def deflection(m, node, wd, tag):
    deck = _model_cards(m) + ["*STEP", "*STATIC", "*CLOAD", f"{node},3,-1.0e4",
                              f"*NSET,NSET=P", f"{node}", "*NODE PRINT,NSET=P", "U", "*END STEP", ""]
    (wd / f"{tag}.inp").write_text("\n".join(deck))
    r = run_ccx(tag, wd, threads=4)
    assert r.returncode == 0, (r.stdout + r.stderr)[-600:]
    txt = (wd / f"{tag}.dat").read_text()
    return -float(re.findall(rf"^\s*{node}\s+\S+\s+\S+\s+(\S+)\s*$", txt, re.M)[-1])


def main() -> None:
    s = shell_surface()
    v60 = build_mesh(s, 0.6).area() * T_REF
    lay = json.load(open(ROOT / "results" / "layouts.json"))
    m, met = ribbed_model(Layout("D", [np.array(x) for x in lay["D"]]), s, v60, h=0.6)
    rid = max(m.mesh.ribs, key=m.rib_length)
    mids = [mid for _, mid, _ in m.mesh.ribs[rid]]
    node = mids[len(mids) // 2]
    rows = [f"point load 10 kN on rib {rid} of D (node {node} at "
            f"{np.round(m.mesh.nodes[node], 2).tolist()}), base shell 50 mm", "rib depth [mm]   w [mm]   w / w_bare"]
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        w0 = None
        for d in (1e-3, 0.05, 0.1, 0.209, 0.4, 0.8):
            mm = dataclasses.replace(m, ribs=[Rib(r.rib_id, r.width, d) for r in m.ribs])
            w = deflection(mm, node, wd, f"d{int(d * 1e3)}")
            w0 = w0 or w
            rows.append(f"{d * 1e3:12.0f}   {w * 1e3:7.4f}   {w / w0:7.3f}")
    OUT.mkdir(exist_ok=True)
    (OUT / "rib_coupling_curved.txt").write_text("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
