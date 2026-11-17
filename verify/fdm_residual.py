"""Equilibrium residual of the form-found shell (force density method, ShellBench).

Vertical equilibrium of every free node of the fixed-plan network:
    q_x (2 z_ij - z_i+1,j - z_i-1,j) + q_y (2 z_ij - z_i,j+1 - z_i,j-1) = w dx dy
The relative residual max|lhs - rhs| / (w dx dy) must be below 1e-6. The bicubic surface used for
meshing must also pass through the network nodes (max |spline - z|).

Usage: python verify/fdm_residual.py   (writes verify/out/fdm_residual.txt)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ribshell.form import shell_surface

OUT = Path(__file__).parent / "out"


def main() -> None:
    s = shell_surface()
    v = s.vault
    z = v.z
    p = v.w * v.dx * v.dy
    lhs = (v.q_x * (2 * z[1:-1, 1:-1] - z[1:-1, 2:] - z[1:-1, :-2])
           + v.q_y * (2 * z[1:-1, 1:-1] - z[2:, 1:-1] - z[:-2, 1:-1]))
    res = np.abs(lhs - p).max() / p
    fit = np.abs(s.z(v.x, v.y) - z).max()
    out = [f"FDM network {z.shape[1]-1} x {z.shape[0]-1}, q_x = {v.q_x:.6g}, q_y = {v.q_y:.6g} N/m, "
           f"rise = {z.max():.6f} m",
           f"max relative vertical residual = {res:.3e}  (limit 1e-6)",
           f"max |spline - network z| = {fit:.3e} m",
           "PASS" if res < 1e-6 and fit < 1e-9 else "FAIL"]
    print("\n".join(out))
    OUT.mkdir(exist_ok=True)
    (OUT / "fdm_residual.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
