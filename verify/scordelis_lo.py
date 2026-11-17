"""Scordelis-Lo roof with S8R shells (brief section 6, first check).

The model is ShellBench's (imported, not copied or edited): whole roof, R = 25, L = 50, 40 deg,
t = 0.25, E = 4.32e8, nu = 0, self-weight 90 per unit area, rigid diaphragms at the curved
ends, free straight edges. The target is the vertical displacement 0.3024 at the free-edge
midpoint (MacNeal & Harder 1985). It runs here with the same ccx build as the P2 study.
Pass: within 2 % on the converged mesh.

Usage: python verify/scordelis_lo.py   (writes verify/out/scordelis_lo.txt)
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from cecore.env import setup_ccx
from shellbench.shellfe import benchmark_row, load_benchmarks

OUT = Path(__file__).parent / "out"
NS = (4, 8, 12, 16)


def main() -> None:
    setup_ccx()   # puts the conda ccx and its DLLs on PATH for ShellBench's lookup
    p = load_benchmarks()["scordelis_lo"]
    rows = []
    with tempfile.TemporaryDirectory() as td:
        for n in NS:
            r = benchmark_row("scordelis_lo", n, Path(td), p)
            rows.append(r)
            print(f"n={n:3d}  u={r['u']:.5f}  u/ref={r['normalised']:.4f}", flush=True)
    last = rows[-1]["normalised"]
    out = [f"Scordelis-Lo roof, S8R (ShellBench model, ccx as used in P2), reference {p['reference_displacement']}"]
    out += [f"{2*r['n']} x {2*r['n']} elements: u = {r['u']:.5f}, u/ref = {r['normalised']:.4f}" for r in rows]
    out.append(f"finest mesh error {100*(last-1):+.2f} %  (limit 2 %)")
    out.append("PASS" if abs(last - 1) < 0.02 else "FAIL")
    print("\n".join(out))
    OUT.mkdir(exist_ok=True)
    (OUT / "scordelis_lo.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
