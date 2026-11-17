"""Zoelly spherical-cap LBA verification.

Classical buckling pressure for a complete thin spherical shell under uniform external pressure
(Zoelly 1915, Timoshenko & Gere 1961 eq. 11-12):

    p_cr = 2 E (t/R)^2 / sqrt(3(1 - nu^2))

For a clamped shallow cap the knockdown is small; pinned edges lower it somewhat.
We mesh a hemisphere (half-angle 90°) with pinned equator and uniform pressure, expecting
the first eigenvalue within ~5 % of the Zoelly formula (the mesh and BC slightly lower it).

Usage:  python verify/zoelly_cap.py
"""
from __future__ import annotations

import math
import tempfile
from pathlib import Path

import gmsh
import numpy as np

from ribshell.constants import E, NU
from ribshell.fe import read_factors

OUT = Path(__file__).parent / "out"


def zoelly_pcr(e: float, nu: float, t: float, r: float) -> float:
    """Classical Zoelly critical pressure [Pa]."""
    return 2 * e * (t / r) ** 2 / math.sqrt(3 * (1 - nu**2))


def hemisphere_mesh(r: float, t: float, h: float) -> tuple[dict, dict, set]:
    """Mesh a hemisphere of radius R, return nodes {id: (x,y,z)}, quads {id: 8-tuple}, edge set."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cap")
    occ = gmsh.model.occ
    # hemisphere surface: sphere cut by z=0 plane
    sph = occ.addSphere(0, 0, 0, r)
    box = occ.addBox(-2 * r, -2 * r, -2 * r, 4 * r, 4 * r, 2 * r)
    occ.cut([(3, sph)], [(3, box)])
    occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.3 * h)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.model.mesh.generate(2)

    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    nodes = {int(t): tuple(p) for t, p in zip(tags, xyz)}

    quads = {}
    for _, tag in gmsh.model.getEntities(2):
        types, etags, enodes = gmsh.model.mesh.getElements(2, tag)
        for ty, et, en in zip(types, etags, enodes):
            if ty != 16:
                continue
            for e, conn in zip(et, en.reshape(-1, 8)):
                c = tuple(int(n) for n in conn)
                # ensure outward normal (radial): signed area in tangent plane at centroid
                ctr = np.mean([nodes[n] for n in c[:4]], axis=0)
                ctr /= np.linalg.norm(ctr)  # unit radial
                p = np.array([nodes[n] for n in c[:4]])
                # project onto tangent plane
                p2d = p - np.outer(p @ ctr, ctr)
                # pick two orthogonal axes in tangent plane
                u = p2d[1] - p2d[0]
                u /= np.linalg.norm(u)
                v = np.cross(ctr, u)
                px = p2d @ u
                py = p2d @ v
                signed = 0.5 * sum(px[i] * py[(i + 1) % 4] - px[(i + 1) % 4] * py[i] for i in range(4))
                if signed < 0:
                    c = (c[0], c[3], c[2], c[1], c[7], c[6], c[5], c[4])
                quads[int(e)] = c
    gmsh.finalize()

    used = {n for c in quads.values() for n in c}
    nodes = {n: p for n, p in nodes.items() if n in used}
    tol = 1e-4 * r
    edge = {n for n, (x, y, z) in nodes.items() if abs(z) < tol}
    return nodes, quads, edge


def cap_deck(nodes, quads, edge, t: float, e: float, nu: float, rho: float, p: float) -> str:
    """CalculiX deck: pinned equator, uniform pressure p (positive = inward on outward-normal shell)."""
    lines = ["*NODE"] + [f"{n},{x:.10g},{y:.10g},{z:.10g}" for n, (x, y, z) in nodes.items()]
    lines += ["*ELEMENT,TYPE=S8R,ELSET=SHELL"]
    lines += [f"{eid}," + ",".join(map(str, c)) for eid, c in quads.items()]
    lines += ["*MATERIAL,NAME=C", "*ELASTIC", f"{e},{nu}", "*DENSITY", f"{rho}"]
    lines += ["*SHELL SECTION,ELSET=SHELL,MATERIAL=C", f"{t}"]
    lines += ["*NSET,NSET=EDGE"] + [",".join(map(str, sorted(edge)[k:k + 16])) for k in range(0, len(edge), 16)]
    lines += ["*BOUNDARY", "EDGE,1,3"]  # pinned
    lines += ["*STEP", "*BUCKLE", "10"]
    lines += ["*DLOAD", f"SHELL,P,{p:.10g}"]
    lines += ["*NODE FILE", "U", "*END STEP", ""]
    return "\n".join(lines)


def main() -> None:
    import os
    import subprocess
    from cecore.env import setup_ccx

    OUT.mkdir(exist_ok=True)
    R, t = 5.0, 0.05  # 5 m radius, 50 mm thick
    h = 0.2  # mesh size (buckle wavelength 2.44*sqrt(Rt) = 1.2 m)
    p_ref = -1e6  # 1 MPa external (inward): normals point outward, +P acts along +normal

    nodes, quads, edge = hemisphere_mesh(R, t, h)
    deck = cap_deck(nodes, quads, edge, t, E, NU, 2500.0, p_ref)

    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "cap.inp").write_text(deck)
        ccx = setup_ccx()
        env = dict(os.environ, OMP_NUM_THREADS="4")
        r = subprocess.run([ccx, "-i", "cap"], cwd=d, env=env, capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(r.stdout + r.stderr)
        lam = read_factors(d / "cap.dat")
        # copy frd for rendering
        (OUT / "cap.frd").write_bytes((d / "cap.frd").read_bytes())

    p_cr_fe = lam[0] * abs(p_ref)
    p_cr_z = zoelly_pcr(E, NU, t, R)
    err = (p_cr_fe - p_cr_z) / p_cr_z * 100

    report = [
        f"Zoelly cap verification (R={R} m, t={t*1e3:.0f} mm, E={E/1e9:.0f} GPa, nu={NU})",
        f"Mesh: {len(quads)} S8R elements, {len(nodes)} nodes, h={h} m",
        f"Zoelly p_cr = {p_cr_z/1e6:.4f} MPa",
        f"FE lambda_1 = {lam[0]:.4f}  =>  p_cr = {p_cr_fe/1e6:.4f} MPa",
        f"Error = {err:+.2f} %",
        "PASS" if abs(err) < 5 else "FAIL",
    ]
    print("\n".join(report))
    (OUT / "zoelly_cap.txt").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
