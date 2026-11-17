#!/usr/bin/env bash
# Full study chain at the study mesh (h = 0.4 m, from verify/mesh_convergence.py). Resumable:
# each step reuses its result file if present. Logs go to results/logs/.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python}
H=${H:-0.4}
export PYTHONPATH=src
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a results/logs/run_all.log; }
step() { local name=$1; shift; log "start $name"; "$PY" "$@" > "results/logs/$name.log" 2>&1; log "done  $name"; }
[ -f "results/variants_h${H}_r1.json" ] || step s02_r1 study/s02_variants.py --h "$H" --vratio 1.0
[ -f "results/E_best_d4_h${H}.json" ]   || step s03_d4 study/s03_optimise_E.py --stage d4 --h "$H"
[ -f "results/E_best_free_h${H}.json" ] || step s03_free study/s03_optimise_E.py --stage free --h "$H"
[ -f "results/gnia_h${H}.json" ]        || step s04 study/s04_gnia.py --h "$H"
step s05 study/s05_volume_sweep.py --h "$H"
for v in ribbed_strip check_fields sway_probe time_buckle; do step "verify_$v" "verify/$v.py"; done
log "all done"
