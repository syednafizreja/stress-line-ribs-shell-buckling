"""Check the stress-resultant extraction (ribshell.fields) on a sagging plate.

Simply supported square plate a = 4 m, t = 100 mm, self-weight q = rho g t, nu = 0.2.
Navier / Timoshenko & Woinowsky-Krieger (1959, Table 8): M_x = M_y = (1 + nu) 0.0368 q a^2 at the
centre (0.0479 q a^2 for nu = 0.3, Table 8). Expect m11 = m22 < 0 (top face in compression under the
'm > 0 puts the top in tension' convention), n ~ 0, and the integration points 5-8 on the top face.

Usage:  python verify/check_fields.py   (writes verify/out/check_fields.txt)
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

from ribshell.constants import G, RHO
from ribshell.fe import Model
from ribshell.fields import BOT_IPS, TOP_IPS, read_ip_stresses, resultants, static_deck
from ribshell.form import Surface
from ribshell.mesh import build_mesh

OUT = Path(__file__).parent / "out"


def flat(a: float) -> Surface:
    g = np.linspace(0, a, 5)
    return Surface(None, RectBivariateSpline(g, g, np.zeros((5, 5)), kx=3, ky=3), a)


def main() -> None:
    from cecore.env import run_ccx

    a, t = 4.0, 0.10
    s = flat(a)
    m = Model(build_mesh(s, 0.1), t)
    cases = {"G": [(G, (0.0, 0.0, -1.0))]}
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "plate.inp").write_text(static_deck(m, cases))
        r = run_ccx("plate", d, threads=4)
        assert r.returncode == 0, r.stdout[-500:]
        ip = read_ip_stresses(d / "plate.dat", 1)[0]
    res = resultants(m, ip, s)
    k = int(np.argmin(np.linalg.norm(res.centroid[:, :2] - a / 2, axis=1)))
    e = int(res.elems[k])
    sxx = ip[e][:, 0]
    q = RHO * G * t
    m_ref = (1 + 0.2) * 0.0368 * q * a**2
    lines = [
        f"plate a={a} m t={t} m h=0.1 m, {len(m.mesh.quads)} S8R; centre element {e} at "
        f"({res.centroid[k,0]:.3f},{res.centroid[k,1]:.3f})",
        f"sxx at integration points 1-8 [kPa]: " + " ".join(f"{v/1e3:8.1f}" for v in sxx),
        f"mean sxx ips {[i+1 for i in TOP_IPS]} = {sxx[list(TOP_IPS)].mean()/1e3:.1f} kPa (expected < 0, top)",
        f"mean sxx ips {[i+1 for i in BOT_IPS]} = {sxx[list(BOT_IPS)].mean()/1e3:.1f} kPa (expected > 0, bottom)",
        f"m11 = {res.m[k,0,0]:.1f}, m22 = {res.m[k,1,1]:.1f}, m12 = {res.m[k,0,1]:.1f} N m/m;"
        f"  reference -{m_ref:.1f}  ratio {-res.m[k,0,0]/m_ref:.4f}",
        f"|n|max = {np.abs(res.n[k]).max():.2f} N/m (expected ~0)",
    ]
    ok = sxx[list(TOP_IPS)].mean() < 0 < sxx[list(BOT_IPS)].mean() and abs(-res.m[k, 0, 0] / m_ref - 1) < 0.03
    lines.append("PASS" if ok else "FAIL")
    print("\n".join(lines))
    OUT.mkdir(exist_ok=True)
    (OUT / "check_fields.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
