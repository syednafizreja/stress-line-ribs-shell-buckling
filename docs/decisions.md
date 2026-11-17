# Methods log: decisions and additional checks

This file lists two kinds of entries. Some are checks added beyond the project brief. Others are modelling choices the author still has to confirm. Each entry names where it lives in the code and what it showed. Entries marked *author to confirm* were set by the analyst, from primary sources where a source exists, and are waiting for the author's sign-off.

## Checks added beyond the brief

| Check | Where | Purpose | Result |
| --- | --- | --- | --- |
| Ribbed slab checked two ways | `verify/ribbed_strip.py` | Compares the study rib model (S8R + B32R with OFFSET2) with an independent model that uses thickened S8R strips. | Beam 0.976, strips 0.999 of the T-beam reference; beam/strip 0.977. To be re-run on the gated mesher. |
| Layout robustness to the field mesh | `verify/layout_robustness.py` | Traces layouts C, D and Dx from fields computed at h = 0.6, 0.4, 0.3 and 0.2 m and measures the Hausdorff distance to the study layout. | C and Dx are stable. D changes topology at h = 0.3 (8 lines instead of 12). This is reported as found, and the 8-line diagonal layout is studied as its own variant, Dd. |
| Mesh convergence with Richardson extrapolation | `verify/mesh_convergence.py` | Chooses the study mesh from the observed order of convergence rather than from a fixed rule. | Pending. |
| Mesh quality gate | `src/ribshell/mesh.py` | Rejects any mesh with min sin(corner angle) < 0.2. This keeps degenerate quads, which give a nonpositive Jacobian, out of the study. | Active for all runs. Earlier results are in `results/superseded/`. |
| Limit-point versus solver-stop classification | `study/s04_gnia.py` | CalculiX has no arc-length method, so each GNIA stop is classified by its end secant stiffness. A solver stop is never reported as a buckling load. | Applied to every GNIA row. |
| Both imperfection signs | `study/s04_gnia.py` | Runs the mode-1 imperfection in both signs. The weaker sign governs. | Applied to every GNIA row. |
| Tool environment fix | `cecore/src/cecore/env.py`, `cecore/tests/test_env.py` | `setup_ccx` grew PATH on every call and crashed long runs (0xC0000017). It now adds only the missing entries, and a test checks this. | Fixed; the test passes. |

## Choices to be confirmed by the author

| Choice | Where | Value and basis |
| --- | --- | --- |
| Stress-line tracing settings | `src/ribshell/layouts.py` | Anisotropy cut-off 0.02 for C. For D and Dx, lines stop below 15 % of the peak moment, where the principal directions are numerical noise. C and D are traced on the D4 group-averaged field with whole orbits. *Author to confirm.* |
| Edge clearance of ribs | `src/ribshell/layouts.py` | Lines whose median distance to the pinned boundary is below 0.5 m are dropped: the edge is already restrained, so such a rib only uses volume. *Author to confirm.* |
| Reference imperfection amplitude | `study/s04_gnia.py` | L/500 (L = 12 m plan span) as the reference, with L/300 and L/1000 as the bounds. To be checked against the code basis adopted for the report. *Author to confirm.* |
| Two-stage optimisation of variant E | `study/s03_optimise_E.py` | Stage 1: D4-symmetric, 6 variables, load directions 0° and 45°. Stage 2: free, 25 variables, 8 load directions, started from the stage-1 optimum. *Author to confirm.* |
| Study mesh size | `verify/mesh_convergence.py` | Set from the convergence study. The actual element edge is about h/2 because of quad subdivision, and it is reported that way. *Author to confirm.* |
| Vendoring `cecore` | `STATUS.md` | Before either repository goes public, vendor `cecore` into it or publish it as a package, so that a fresh clone installs. *Author to confirm.* |
