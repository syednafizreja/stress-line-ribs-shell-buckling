"""Variant E: CMA-ES search for the zonal thickness that maximises the worst-direction lambda_c.

Objective (maximised): min over load directions theta of lambda_c(theta), G + 0.2 g lateral,
lowest positive factor of 10 modes (so a switch of the critical mode is handled by the min).
Volume is held at V60 exactly by the projection in ribshell.zonal; CMA-ES searches in
x = log(u), so the design is positive and scale-free.

Two stages (see docs/decisions.md; author to confirm):
  1. "d4":   6 variables (D4 orbits of a 5 x 5 zoning), theta in {0, 45}: by symmetry these two
             directions give the minimum over all eight (checked in s02 for every D4 layout).
  2. "free": 25 variables, no symmetry, all eight directions, warm-started at the stage-1 optimum.
             It tests whether breaking symmetry helps under the four-direction sway.

Every evaluation is appended to results/E_<stage>_h<h>.jsonl (the run can be resumed: logged
designs are not recomputed). Best design -> results/E_best_<stage>_h<h>.json.

Usage: python study/s03_optimise_E.py --stage d4 --h 0.6
       python study/s03_optimise_E.py --stage free --h 0.6
"""
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import cma
import numpy as np

from ribshell.constants import T_REF
from ribshell.fe import Model, run_buckle
from ribshell.form import shell_surface
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
NZ = 5
WORKERS, THREADS = 4, 4


def lowest_positive(f):
    pos = [x for x in f if x > 0]
    return min(pos) if pos else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("d4", "free"), required=True)
    ap.add_argument("--h", type=float, default=0.6)
    ap.add_argument("--maxfevals", type=int, default=None)
    a = ap.parse_args()
    s = shell_surface()
    mesh = build_mesh(s, a.h)
    v60 = mesh.area() * T_REF
    des = ZonalDesign(mesh, NZ, v60, symmetry="d4" if a.stage == "d4" else None)
    angles = (0, 45) if a.stage == "d4" else tuple(range(0, 360, 45))
    tag = f"{a.stage}_h{a.h:g}"
    log_path = RES / f"E_{tag}.jsonl"
    cache = {}
    if log_path.exists():
        for ln in log_path.read_text().splitlines():
            r = json.loads(ln)
            cache[tuple(np.round(r["x"], 10))] = r

    if a.stage == "d4":
        x0, sigma0, popsize, maxfe = np.zeros(des.n_vars), 0.3, 10, a.maxfevals or 300
    else:
        d4 = json.load(open(RES / f"E_best_d4_h{a.h:g}.json"))
        x0 = np.log(np.asarray(d4["zone_thickness_m"]))   # a valid u (the projection is scale-free)
        sigma0, popsize, maxfe = 0.15, 14, a.maxfevals or 560

    wd = Path(tempfile.mkdtemp(prefix=f"E_{tag}_"))   # run_buckle removes each run's files
    counter = [len(cache)]

    def lam(job):
        k, ang, tz = job
        m = Model(mesh, des.element_thickness(tz), surface=s)
        return lowest_positive(run_buckle(m, wd, f"e{k}_{ang}", "combined", ang, threads=THREADS).factors)

    def evaluate_batch(X):
        out, todo = [None] * len(X), []
        for i, x in enumerate(X):
            key = tuple(np.round(x, 10))
            if key in cache:
                out[i] = cache[key]
            else:
                todo.append(i)
        tzs = {i: des.zone_thickness(np.exp(X[i])) for i in todo}
        jobs = []
        for i in todo:
            counter[0] += 1
            jobs += [(counter[0] * 100 + i, ang, tzs[i]) for ang in angles]
        with ThreadPoolExecutor(WORKERS) as ex:
            vals = list(ex.map(lam, jobs))
        it = iter(vals)
        with open(log_path, "a") as fh:
            for i in todo:
                lc = {ang: next(it) for ang in angles}
                r = {"x": list(map(float, X[i])), "zone_thickness_m": tzs[i].tolist(),
                     "lambda_c_by_angle": lc, "objective": min(lc.values())}
                cache[tuple(np.round(X[i], 10))] = r
                fh.write(json.dumps(r) + "\n")
                out[i] = r
        return out

    # reference: the uniform design (x = 0 for d4) or the warm start (free)
    ref = evaluate_batch([x0])[0]
    print(f"start objective {ref['objective']:.3f}", flush=True)
    es = cma.CMAEvolutionStrategy(x0, sigma0, {"popsize": popsize, "maxfevals": maxfe, "seed": 20261117,
                                               "verbose": -9, "tolfun": 1e-3, "tolx": 1e-4})
    best = ref
    while not es.stop():
        X = es.ask()
        R = evaluate_batch(X)
        es.tell(X, [-r["objective"] for r in R])
        gen_best = max(R, key=lambda r: r["objective"])
        if gen_best["objective"] > best["objective"]:
            best = gen_best
        print(f"gen {es.countiter:3d}  fevals {es.countevals:4d}  gen-best {gen_best['objective']:.3f}"
              f"  best {best['objective']:.3f}  sigma {es.sigma:.3f}", flush=True)
    best = {**best, "stage": a.stage, "h": a.h, "nz": NZ, "V_m3": v60,
            "volume_check_m3": float(des.zone_area @ np.asarray(best["zone_thickness_m"])),
            "uniform_or_start_objective": ref["objective"], "stop": {k: str(v) for k, v in es.stop().items()},
            "fevals": es.countevals}
    json.dump(best, open(RES / f"E_best_{tag}.json", "w"), indent=1)
    shutil.rmtree(wd, ignore_errors=True)
    print(json.dumps({k: best[k] for k in ("objective", "uniform_or_start_objective", "fevals")}))
    print("zone t [mm]:\n", (np.asarray(best["zone_thickness_m"]).reshape(NZ, NZ) * 1e3).round(1))


if __name__ == "__main__":
    main()
