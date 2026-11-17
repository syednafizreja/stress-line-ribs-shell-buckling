"""Volume sweep V / V60 = 1.0, 1.1, ..., 1.5 for every variant, with A1-A3 carbon.

For each volume ratio r:
  * A, B, C, D, Dd, Dx: s02 at V = r x V60 (results/variants_h<h>_r<r>.json);
  * E: the D4-symmetric CMA-ES stage of s03 at the same V (results/E_best_d4_h<h>[_r<r>].json).
    The free 25-zone stage is run at r = 1 only (s03 --stage free), as a check that symmetry
    does not cap E; its gain over the D4 optimum is reported alongside.
Existing result files are reused, so the sweep can be resumed.

Carbon: A1-A3 only, one factor for in-situ concrete. IStructE (2022) "How to calculate embodied
carbon", 2nd ed., Table 2.3, UK in-situ C32/40, 25 % GGBS: default 0.120 kgCO2e/kg (typical
bounds 0.063 to 0.149), from the ICE database (Circular Ecology 2020). The table has no C30/37
row; C32/40 is the next class up, so the value is on the high side for C30/37. Reinforcement is
excluded (scope section 3). Equal volume means equal carbon, so carbon is a relabelled volume axis.

Recovery (the hypothesis test): R_X = (lambda_X - lambda_A) / (lambda_E - lambda_A).

Outputs: results/sweep_h<h>.json, results/sweep_h<h>.txt
Usage  : python study/s05_volume_sweep.py [--h 0.6]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from s02_variants import run as run_variants  # noqa: E402

from ribshell.constants import RHO  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RATIOS = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
ECF = {"default": 0.120, "low": 0.063, "high": 0.149}   # kgCO2e per kg, A1-A3
RIBBED = ("B", "C", "D", "Dd", "Dx")


def e_file(h: float, r: float, stage: str = "d4") -> Path:
    return RES / f"E_best_{stage}_h{h:g}{'' if r == 1 else f'_r{r:g}'}.json"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.6)
    h = ap.parse_args().h
    rows = []
    for r in RATIOS:
        vf = RES / f"variants_h{h:g}_r{r:g}.json"
        var = json.load(open(vf)) if vf.exists() else run_variants(h, r)
        if not e_file(h, r).exists():
            subprocess.run([sys.executable, str(ROOT / "study" / "s03_optimise_E.py"), "--stage", "d4",
                            "--h", str(h), "--vratio", str(r)], check=True)
        e = json.load(open(e_file(h, r)))
        lam = {k: v["lambda_c"] for k, v in var["variants"].items()}
        lam["E"] = e["objective"]
        mass = var["V_m3"] * RHO
        row = {"vratio": r, "V_m3": var["V_m3"], "mean_t_mm": var["V_m3"] / (var["V60_m3"] / 0.060) * 1e3,
               "tCO2e": {k: mass * f / 1e3 for k, f in ECF.items()}, "lambda_c": lam,
               "lambda_g": {**{k: v["lambda_g"] for k, v in var["variants"].items()}},
               "recovery": {k: (lam[k] - lam["A"]) / (lam["E"] - lam["A"]) for k in RIBBED if k in lam},
               "fabrication": {k: {f: var["variants"][k][f] for f in ("rib_length_m", "crossings", "rib_depth_mm",
                                                                       "min_radius_m", "n_ribs")}
                               for k in RIBBED if k in var["variants"]},
               "E_thickness_mm": [min(e["zone_thickness_m"]) * 1e3, max(e["zone_thickness_m"]) * 1e3]}
        ef = e_file(h, r, "free")
        if ef.exists():
            row["E_free"] = json.load(open(ef))["objective"]
        rows.append(row)
        print(r, {k: round(v, 1) for k, v in lam.items()}, flush=True)
    json.dump({"h": h, "ecf_kgco2e_per_kg": ECF, "ecf_source": "IStructE 2022, HTCEC 2nd ed., Table 2.3 (ICE 2020)",
               "rows": rows}, open(RES / f"sweep_h{h:g}.json", "w"), indent=1)
    names = ["A", *[k for k in RIBBED if k in rows[0]["lambda_c"]], "E"]
    out = [f"Volume sweep, h = {h} m, lambda_c = min over 8 lateral directions (E: d4 optimum), a_h = 0.2 g",
           "r     V[m3]  tCO2e  " + "".join(f"{k:>8s}" for k in names) + "   " +
           "".join(f"  R_{k:<4s}" for k in names[1:-1])]
    for row in rows:
        out.append(f"{row['vratio']:.1f}  {row['V_m3']:6.2f}  {row['tCO2e']['default']:5.2f}  " +
                   "".join(f"{row['lambda_c'][k]:8.1f}" for k in names) + "   " +
                   "".join(f"  {row['recovery'][k]:6.3f}" for k in names[1:-1]))
    out.append(f"carbon factor {ECF['default']} kgCO2e/kg (band {ECF['low']}-{ECF['high']}), A1-A3, no reinforcement")
    (RES / f"sweep_h{h:g}.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
