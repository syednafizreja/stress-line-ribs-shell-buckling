"""Layout Dd: the diagonal-only moment-line layout.

The seismic moment field traced on the h = 0.3 m field mesh gives 8 diagonal lines instead of the
12 lines found at h = 0.6, 0.4 and 0.2 m (verify/layout_robustness.py). Rather than pick one, both
are studied: D (12 lines, the h = 0.4 trace) and Dd (8 diagonal lines, the h = 0.3 trace).

Input : verify/out/layout_robustness.json      Output: "Dd" added to results/layouts.json
Usage : python study/s01b_layout_Dd.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    rob = json.load(open(ROOT / "verify" / "out" / "layout_robustness.json"))
    dd = rob["D"]["0.3"]["lines"]
    lay_path = ROOT / "results" / "layouts.json"
    lay = json.load(open(lay_path))
    lay["Dd"] = dd
    json.dump(lay, open(lay_path, "w"))
    print(f"Dd: {len(dd)} lines added to {lay_path.name}")


if __name__ == "__main__":
    main()
