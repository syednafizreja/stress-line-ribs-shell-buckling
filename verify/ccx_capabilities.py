"""Smoke tests of the CalculiX features Project 2 depends on (run before any study code).

T1  *BUCKLE after a *STATIC preload: what does the buckling factor scale?
    Pinned-pinned B32R column, closed form P_cr = pi^2 EI / L^2.
T2  B32R rib with OFFSET sharing nodes with S8R shell: does it run, and is the
    midspan deflection of a simply supported T-strip near the hand T-beam value?
T3  Arc-length (Riks) for the GNIA limit point: is it accepted by this build?

Findings (17 Nov 2026, ccx 2.23; manual = CalculiX 2.23 User's Manual, G. Dhondt, sec. 7.5, 7.3, 6.2):
- *BUCKLE drops all earlier loads unless the buckle step is *STEP, PERTURBATION (manual 7.5).
  With PERTURBATION the last *STATIC load is a fixed preload and lambda scales only the
  buckle-step load (variants d, e). Use: *STATIC G, then *STEP, PERTURBATION / *BUCKLE with a_h.
- OFFSET2 > 0 puts the beam axis on the -n2 side (manual 7.3: OFFSET2=0.5 makes the node line the
  top surface). With n2 = +z, a positive OFFSET2 hangs the rib below the shell (frd z range).
- S8R normal: nodes clockwise when looking along the normal (manual 6.2); a positive P load acts
  along +normal in this build (T2, P+ gives +uz). Gravity-like pressure on a +z-normal shell is P < 0.
- The manual has no Riks or arc-length entry; *STATIC, RIKS is ignored with a warning (T3).

Usage:  python verify/ccx_capabilities.py      (writes verify/out/ccx_capabilities.txt)
"""
from __future__ import annotations

import math
import re
import tempfile
from pathlib import Path

from cecore.env import ccx_version, run_ccx, setup_ccx

E, NU = 33e9, 0.2
OUT = Path(__file__).parent / "out"


# ---------------------------------------------------------------- T1 column
def column_deck(pre: float | None, buck: float, perturbation: bool = False,
                nlgeom: bool = False) -> str:
    L, b, n_el = 3.0, 0.1, 10
    n_nodes = 2 * n_el + 1
    s = ["*NODE"] + [f"{i+1}, 0, 0, {L*i/(n_nodes-1):.6f}" for i in range(n_nodes)]
    s += ["*ELEMENT, TYPE=B32R, ELSET=COL"]
    s += [f"{e+1}, {2*e+1}, {2*e+2}, {2*e+3}" for e in range(n_el)]
    s += ["*MATERIAL, NAME=C", "*ELASTIC", f"{E}, {NU}",
          "*BEAM SECTION, ELSET=COL, MATERIAL=C, SECTION=RECT", f"{b}, {b}", "1, 0, 0",
          "*BOUNDARY", "1, 1, 3", "1, 6, 6", f"{n_nodes}, 1, 2"]
    if pre is not None:
        s += ["*STEP, NLGEOM" if nlgeom else "*STEP", "*STATIC", "*CLOAD", f"{n_nodes}, 3, {-pre}",
              "*END STEP"]
    s += ["*STEP, PERTURBATION" if perturbation else "*STEP", "*BUCKLE", "3", "*CLOAD", f"{n_nodes}, 3, {-buck}", "*END STEP"]
    return "\n".join(s) + "\n"


def buckling_factors(dat: str) -> list[float]:
    blk = dat.split("B U C K L I N G")[-1]
    return [float(m.group(2)) for m in re.finditer(r"^\s+(\d+)\s+([-+0-9.Ee]+)\s*$", blk, re.M)]


# ---------------------------------------------------------------- T2 T-strip
def tstrip_deck(nx: int = 24, ny: int = 4, offset_sign: float = 1.0,
                p: float = 10e3) -> tuple[str, dict]:
    """Simply supported strip, span 3 m, flange 0.6 m x 50 mm, rib 150 x 150 mm below.

    S8R flange (8-node quads), B32R rib along the centre line y = 0, offset so the
    rib hangs below the shell mid-surface. Uniform pressure on the flange.
    """
    Ls, W, tf, bw, hw = 3.0, 0.6, 0.05, 0.15, 0.15
    # node grid with mid-side nodes: (2nx+1) x (2ny+1)
    NX, NY = 2 * nx + 1, 2 * ny + 1
    nid = lambda i, j: j * NX + i + 1  # noqa: E731
    nodes = [f"{nid(i,j)}, {Ls*i/(NX-1):.6f}, {-W/2 + W*j/(NY-1):.6f}, 0"
             for j in range(NY) for i in range(NX)
             if not (i % 2 == 1 and j % 2 == 1)]  # S8R has no centre node
    els = []
    for jj in range(ny):
        for ii in range(nx):
            i0, j0 = 2 * ii, 2 * jj
            c = [nid(i0, j0), nid(i0+2, j0), nid(i0+2, j0+2), nid(i0, j0+2),
                 nid(i0+1, j0), nid(i0+2, j0+1), nid(i0+1, j0+2), nid(i0, j0+1)]
            els.append(f"{len(els)+1}, " + ", ".join(map(str, c)))
    jc = ny  # centre row index in node grid
    beams = [f"{10000+e+1}, {nid(2*e, jc)}, {nid(2*e+1, jc)}, {nid(2*e+2, jc)}" for e in range(nx)]
    # offset of the rib centroid below shell mid-surface: tf/2 + hw/2, in units of rib height
    off = (tf / 2 + hw / 2) / hw
    s = ["*NODE"] + nodes
    s += ["*ELEMENT, TYPE=S8R, ELSET=FL"] + els
    s += ["*ELEMENT, TYPE=B32R, ELSET=RIB"] + beams
    s += ["*MATERIAL, NAME=C", "*ELASTIC", f"{E}, 0.0",
          "*SHELL SECTION, ELSET=FL, MATERIAL=C", f"{tf}",
          f"*BEAM SECTION, ELSET=RIB, MATERIAL=C, SECTION=RECT, OFFSET2={offset_sign*off:.6f}",
          f"{bw}, {hw}", "0, 1, 0"]
    # supports: line x=0 and x=Ls, vertical; one point fixed in x, y
    s += ["*NSET, NSET=SUP"] + [f"{nid(0,j)}, {nid(NX-1,j)}," for j in range(NY) if j % 2 == 0]
    s += ["*BOUNDARY", "SUP, 3, 3", f"{nid(0,jc)}, 1, 2", f"{nid(NX-1,jc)}, 2, 2"]
    s += ["*STEP", "*STATIC", "*DLOAD"] + [f"{e+1}, P, {p}" for e in range(len(els))]
    s += ["*NODE PRINT, NSET=MID", "U", "*NODE FILE", "U", "*END STEP"]
    s.insert(s.index("*BOUNDARY"), "*NSET, NSET=MID")
    s.insert(s.index("*BOUNDARY"), f"{nid(nx, jc)}")
    # hand T-beam (E with nu=0, full flange width effective for a 0.6 m flange on 3 m span)
    Af, Aw = W * tf, bw * hw
    zf, zw = 0.0, -(tf / 2 + hw / 2)
    zc = (Af * zf + Aw * zw) / (Af + Aw)
    I = W * tf**3 / 12 + Af * (zf - zc)**2 + bw * hw**3 / 12 + Aw * (zw - zc)**2
    q = abs(p) * W
    w_hand = 5 * q * Ls**4 / (384 * E * I)
    return "\n".join(s) + "\n", {"w_hand": w_hand, "mid": nid(nx, jc)}


def midspan_uz(dat: str, node: int) -> float:
    m = re.search(rf"^\s*{node}\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", dat, re.M)
    return float(m.group(3)) if m else float("nan")


def frd_z_range(frd: str) -> tuple[float, float]:
    """z range of all nodes in the .frd, which holds the expanded 3-D nodes of shells and beams."""
    zs = [float(ln[37:49]) for ln in frd.splitlines() if ln.startswith(" -1") and len(ln) >= 49]
    return (min(zs), max(zs)) if zs else (float("nan"), float("nan"))


# ---------------------------------------------------------------- T3 Riks
def riks_deck() -> str:
    base = column_deck(None, 1.0).split("*STEP")[0]
    return base + "*STEP, NLGEOM\n*STATIC, RIKS\n*CLOAD\n21, 3, -1e5\n*END STEP\n"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    lines = [f"ccx: {setup_ccx()}", f"version: {ccx_version()}"]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)

        def run(name: str, deck: str) -> tuple[int, str, str]:
            (d / f"{name}.inp").write_text(deck)
            r = run_ccx(name, d, timeout=600)
            dat = (d / f"{name}.dat").read_text() if (d / f"{name}.dat").exists() else ""
            frd = d / f"{name}.frd"
            run.frd = frd.read_text() if frd.exists() else ""
            return r.returncode, dat, r.stdout + r.stderr

        # T1
        EI = E * 0.1**4 / 12
        pcr = math.pi**2 * EI / 3.0**2
        lines.append(f"\nT1 column: closed-form P_cr = {pcr/1e3:.2f} kN")
        res = {}
        for tag, pre, buck in (("a_nopre", None, 100e3), ("b_pre100_buck100", 100e3, 100e3),
                               ("c_pre100_buck50", 100e3, 50e3)):
            rc, dat, _ = run(f"t1{tag}", column_deck(pre, buck))
            lam = buckling_factors(dat)
            res[tag] = lam[0] if lam else float("nan")
            lines.append(f"  {tag}: rc={rc} lambda1={res[tag]:.4f}  lambda*P_buck={res[tag]*buck/1e3:.2f} kN")
        a, b, c = res["a_nopre"], res["b_pre100_buck100"], res["c_pre100_buck50"]
        lines.append(f"  no-preload error vs closed form: {100*(a*100e3/pcr-1):+.2f} %")
        if abs(c * 50e3 + 100e3 - pcr) / pcr < 0.03:
            lines.append("  RESULT: preload is kept in the stiffness; lambda scales only the buckle-step load "
                         "-> put G in the static step and a_h in the buckle step, or G+a_h both in the buckle step.")
        elif abs(c * 50e3 - pcr) / pcr < 0.03:
            lines.append("  RESULT (a-c): without PERTURBATION the preload is dropped; lambda scales the buckle-step "
                         "load alone.")
        else:
            lines.append("  RESULT: unclear -> inspect decks by hand.")

        # T1 extra: PERTURBATION on the buckle step (manual 7.5: preload kept only then),
        # with a linear and an NLGEOM preload step. Kept preload predicts lambda*50 + 100 = P_cr.
        for tag, nl in (("d_pert_pre100_buck50", False), ("e_pert_nlgeom_pre100_buck50", True)):
            rc, dat, _ = run(f"t1{tag}", column_deck(100e3, 50e3, perturbation=True, nlgeom=nl))
            lam = buckling_factors(dat)
            v = lam[0] if lam else float("nan")
            lines.append(f"  {tag}: rc={rc} lambda1={v:.4f}  lambda*50+100 = {v*50+100:.2f} kN "
                         f"(error vs P_cr {100*((v*50e3+100e3)/pcr-1):+.2f} %)")

        # T2
        lines.append("\nT2 S8R + offset B32R T-strip (S8R nodes counter-clockwise from above -> normal +z;"
                     " beam along +x, n1=(0,1,0) -> n2 = t x n1 = +z)")
        for tag, sgn, p in (("offset-, P+", -1.0, 10e3), ("offset+, P+", 1.0, 10e3),
                            ("offset+, P-", 1.0, -10e3)):
            deck, info = tstrip_deck(offset_sign=sgn, p=p)
            rc, dat, log = run("t2strip", deck)
            uz = midspan_uz(dat, info["mid"])
            z0, z1 = frd_z_range(run.frd)
            lines.append(f"  {tag}: rc={rc} uz_mid={uz*1e3:+.4f} mm, |hand| {abs(info['w_hand'])*1e3:.4f} mm, "
                         f"ratio {abs(uz/info['w_hand']):.3f}; frd node z range [{z0:+.3f}, {z1:+.3f}] m")
            if rc != 0:
                lines.append("  ccx log tail: " + " | ".join(log.strip().splitlines()[-4:]))

        # T3
        rc, dat, log = run("t3riks", riks_deck())
        tail = " | ".join([ln for ln in log.splitlines() if "ERROR" in ln.upper() or "WARN" in ln.upper()][:3])
        lines.append(f"\nT3 *STATIC, RIKS: rc={rc}; messages: {tail or '(none)'}")
    txt = "\n".join(lines) + "\n"
    (OUT / "ccx_capabilities.txt").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
