"""Reproduce every result, figure and export of the study from scratch.

    conda env create -f environment.yml && conda activate research
    python reproduce.py            # full chain at the study mesh h = 0.4 m (about 12 h on 16 threads)
    python reproduce.py --quick    # figures and exports only, from the result files in results/

Order: verification checks (verify/), layouts (s01, s01b), variants at V60 (s02), optimised
thickness E (s03, D4 then free stage), nonlinear check (s04), volume sweep (s05), lateral-acceleration sensitivity (s08),
exports and mode shapes (s07), figures (s06). For a resumable run of the study steps use study/run_all.sh.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
H = "0.4"
VERIFY = ["smoke_env", "ccx_capabilities", "scordelis_lo", "zoelly_cap", "ribbed_strip", "fdm_residual",
          "check_fields", "layout_robustness", "mesh_convergence", "rib_effectiveness", "rib_coupling_curved", "sway_probe", "time_buckle"]
STUDY = [["s01_fields_layouts.py"], ["s01b_layout_Dd.py"], ["s02_variants.py", "--h", H, "--vratio", "1.0"],
         ["s03_optimise_E.py", "--stage", "d4", "--h", H], ["s03_optimise_E.py", "--stage", "free", "--h", H],
         ["s04_gnia.py", "--h", H], ["s05_volume_sweep.py", "--h", H],
         ["s08_ah_sensitivity.py", "--h", H]]


def run(args: list[str]) -> None:
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    print(">", " ".join(args), flush=True)
    subprocess.run([sys.executable, *args], cwd=ROOT, env=env, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="only exports and figures from existing results")
    a = ap.parse_args()
    if not a.quick:
        for v in VERIFY:
            if (ROOT / "verify" / f"{v}.py").exists():
                run([f"verify/{v}.py"])
        for s in STUDY:
            if (ROOT / "study" / s[0]).exists():
                run([f"study/{s[0]}", *s[1:]])
    run(["study/s07_exports.py", "--h", H] + (["--skip-modes"] if a.quick else []))
    run(["study/s06_figures.py", "--h", H])


if __name__ == "__main__":
    main()
