"""Sensitivity of lambda_c to the lateral acceleration, at the registered volume.

The designs are the recorded V60 models (h = 0.4 m). Changing a_h changes only the
load. E is the D4 optimum found at 0.20 g; it is not re-optimised. Each new case is
G + a_h together, eight directions, lowest positive factor of 10 modes.

The 0.20 g row is copied from results/variants_h<h>_r1.json and
results/E_best_d4_h<h>.json. Angles 0 and 45 of A are checked against
verify/out/sway_probe.txt, which analysed the same uniform shell.

Writes results/ah_h<h>.json and results/ah_h<h>.txt
Usage:  python study/s08_ah_sensitivity.py [--h 0.4]
"""
from __future__ import annotations

import argparse
import json
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from cecore.env import setup_ccx
from ribshell.constants import T_REF
from ribshell.fe import Model, run_buckle
from ribshell.form import shell_surface
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ANGLES = tuple(range(0, 360, 45))
A_HS = (0.15, 0.30)
WORKERS, THREADS = 3, 4
PROBE = {
    0.15: {0: 348.65, 45: 353.75},
    0.30: {0: 327.04, 45: 336.97},
}


def lowest_positive(factors: list[float]) -> float:
    pos = [x for x in factors if x > 0]
    if not pos:
        raise RuntimeError(f"no positive buckling factor in {factors}")
    return min(pos)


def models(h: float) -> dict[str, Model]:
    s = shell_surface()
    mesh = build_mesh(s, h)
    v60 = mesh.area() * T_REF
    out: dict[str, Model] = {"A": Model(mesh, T_REF, surface=s)}
    lay = json.load(open(RES / "layouts.json"))
    for n in ("B", "C", "D", "Dd", "Dx"):
        if n in lay:
            out[n], _ = ribbed_model(Layout(n, [np.array(x) for x in lay[n]]), s, v60, h=h)
    e = json.load(open(RES / f"E_best_d4_h{h:g}.json"))
    des = ZonalDesign(mesh, e["nz"], v60)
    tz = des.zone_thickness(np.asarray(e["zone_thickness_m"], dtype=float))
    out["E"] = Model(mesh, des.element_thickness(tz), surface=s)
    return out


def recorded(h: float) -> dict[str, dict]:
    var = json.load(open(RES / f"variants_h{h:g}_r1.json"))["variants"]
    e = json.load(open(RES / f"E_best_d4_h{h:g}.json"))
    row = {}
    for n, r in var.items():
        row[n] = {"lambda_c": r["lambda_c"], "worst_angle": r["worst_angle"],
                  "lambda_c_by_angle": {str(k): v for k, v in r["lambda_c_by_angle"].items()},
                  "source": f"variants_h{h:g}_r1.json"}
    row["E"] = {"lambda_c": e["objective"],
                "worst_angle": int(min(e["lambda_c_by_angle"], key=lambda k: e["lambda_c_by_angle"][k])),
                "lambda_c_by_angle": {str(k): v for k, v in e["lambda_c_by_angle"].items()},
                "source": f"E_best_d4_h{h:g}.json"}
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.4)
    h = ap.parse_args().h
    setup_ccx()
    ms = models(h)
    jobs = [(n, a_h, ang) for n in ms for a_h in A_HS for ang in ANGLES]
    wd = Path(tempfile.mkdtemp(prefix="ah_"))

    def one(job: tuple[str, float, int]) -> tuple[str, float, int, float]:
        n, a_h, ang = job
        tag = f"{n}_a{int(round(a_h * 100))}_c{ang}"
        sub = wd / tag
        fac = lowest_positive(run_buckle(ms[n], sub, "job", "combined", ang, a_h=a_h, threads=THREADS).factors)
        print(f"{n} a_h={a_h:.2f} {ang:3d}  {fac:.2f}", flush=True)
        return n, a_h, ang, fac

    with ThreadPoolExecutor(WORKERS) as ex:
        done = list(ex.map(one, jobs))

    by: dict[tuple[str, float], dict[int, float]] = {}
    for n, a_h, ang, fac in done:
        by.setdefault((n, a_h), {})[ang] = fac

    fresh = {}
    for (n, a_h), lc in by.items():
        worst = min(lc, key=lc.get)
        fresh.setdefault(str(a_h), {})[n] = {
            "lambda_c": lc[worst], "worst_angle": worst,
            "lambda_c_by_angle": {str(a): lc[a] for a in ANGLES},
        }

    checks = []
    ok = True
    for a_h, expect in PROBE.items():
        got = fresh[str(a_h)]["A"]["lambda_c_by_angle"]
        for ang, ref in expect.items():
            diff = abs(got[str(ang)] - ref)
            checks.append({"a_h": a_h, "angle": ang, "probe": ref, "s08": got[str(ang)], "abs_diff": diff})
            if diff > 0.05:
                ok = False

    out = {
        "h": h, "V": "V60", "a_h_new": list(A_HS),
        "note": "E is the recorded D4 optimum at 0.20 g, reanalysed, not re-optimised.",
        "recorded_0.20": recorded(h),
        "rows": fresh,
        "probe_check_A": checks,
        "probe_check_ok": ok,
    }
    RES.mkdir(exist_ok=True)
    json.dump(out, open(RES / f"ah_h{h:g}.json", "w"), indent=1)
    names = [n for n in ("A", "B", "C", "D", "Dd", "Dx", "E") if n in ms]
    lines = [
        f"Lateral acceleration sensitivity, h = {h} m, V = V60, lambda_c = min over 8 directions",
        "E is the D4 optimum recorded at 0.20 g. It is reanalysed at the other accelerations, not re-optimised.",
        "a_h     " + " ".join(f"{n:>8s}" for n in names),
    ]
    rec = out["recorded_0.20"]
    lines.append("0.20  " + " ".join(f"{rec[n]['lambda_c']:8.2f}" for n in names) + "   recorded")
    for a_h in A_HS:
        lines.append(f"{a_h:.2f}  " + " ".join(f"{fresh[str(a_h)][n]['lambda_c']:8.2f}" for n in names))
    lines.append("")
    lines.append("Check of A at 0 and 45 deg against verify/out/sway_probe.txt (limit 0.05):")
    for c in checks:
        lines.append(f"  a_h {c['a_h']:.2f}  {c['angle']:3d}  probe {c['probe']:.2f}  s08 {c['s08']:.2f}  |d| {c['abs_diff']:.4f}")
    lines.append("probe check passed" if ok else "probe check FAILED")
    text = "\n".join(lines) + "\n"
    (RES / f"ah_h{h:g}.txt").write_text(text)
    print(text, end="")
    if not ok:
        raise SystemExit(1)
    print("ah sensitivity written", flush=True)


if __name__ == "__main__":
    main()
