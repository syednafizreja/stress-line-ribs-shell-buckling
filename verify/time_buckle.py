"""Time one *BUCKLE run (10 modes) on the form-found 12 x 12 m shell at three mesh sizes.

Records wall time, peak RAM of ccx, solver, element count and lambda_1 for the lateral (G preload,
a_h scaled) and combined definitions. Sets the compute budget for CMA-ES and the volume sweep.

Usage:  python verify/time_buckle.py      (writes verify/out/time_buckle.txt)
"""
from __future__ import annotations

import os
import subprocess
import tempfile
import threading
import time
from pathlib import Path

import psutil
from cecore.env import setup_ccx

from ribshell.constants import T_REF
from ribshell.fe import Model, buckle_deck, read_factors
from ribshell.form import shell_surface
from ribshell.mesh import build_mesh

OUT = Path(__file__).parent / "out"


def run_measured(deck: str, d: Path, name: str, threads: int) -> tuple[float, float, str, list[float]]:
    (d / f"{name}.inp").write_text(deck)
    ccx = setup_ccx()  # sets PATH first (DLLs), so build env after it
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    t0 = time.perf_counter()
    p = subprocess.Popen([ccx, "-i", name], cwd=d, env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True)
    peak = [0.0]

    def watch():
        try:
            pp = psutil.Process(p.pid)
            while p.poll() is None:
                rss = pp.memory_info().rss + sum(c.memory_info().rss for c in pp.children(recursive=True))
                peak[0] = max(peak[0], rss)
                time.sleep(0.2)
        except psutil.Error:
            pass

    th = threading.Thread(target=watch, daemon=True)
    th.start()
    log, _ = p.communicate()
    wall = time.perf_counter() - t0
    th.join(1)
    if p.returncode:
        raise RuntimeError(log[-800:])
    solver = next((s for s in ("PARDISO", "PaStiX", "SPOOLES", "TAUCS") if s.lower() in log.lower()), "?")
    return wall, peak[0] / 2**20, solver, read_factors(d / f"{name}.dat")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    surf = shell_surface()
    lines = ["Buckle timing, uniform t = 60 mm, pinned edges, 10 modes, 4 threads"]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        for h in (0.6, 0.4, 0.3, 0.2):
            mesh = build_mesh(surf, h)
            m = Model(mesh, T_REF)
            for kind in ("lateral", "combined"):
                wall, ram, solver, lam = run_measured(buckle_deck(m, kind, 0.0), d, f"h{int(h*100)}_{kind}", 4)
                lines.append(f"h={h:.2f} m  elems={len(mesh.quads):5d} nodes={len(mesh.nodes):6d}  {kind:8s} "
                             f"wall={wall:6.1f} s  peakRAM={ram:7.0f} MB  solver={solver}  "
                             f"lambda1={lam[0] if lam else float('nan'):.4f}")
                print(lines[-1], flush=True)
    (OUT / "time_buckle.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
