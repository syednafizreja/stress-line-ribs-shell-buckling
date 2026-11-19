# Methods log: decisions and additional checks

This file lists two kinds of entries. Some are checks added beyond the scope document (docs/scope.md). Others are modelling choices the author still has to confirm. Each entry names where it lives in the code and what it showed. Entries marked *author to confirm* were set by the analyst, from primary sources where a source exists, and are waiting for the author's sign-off.

## Checks added beyond the scope document

| Check | Where | Purpose | Result |
| --- | --- | --- | --- |
| Ribbed slab checked two ways | `verify/ribbed_strip.py` | Compares the study rib model (S8R + B32R with OFFSET2) with an independent model that uses thickened S8R strips. | Beam 0.976, strips 0.999 of the T-beam reference; beam/strip 0.977. To be re-run on the gated mesher. |
| Layout robustness to the field mesh | `verify/layout_robustness.py` | Traces layouts C, D and Dx from fields computed at h = 0.6, 0.4, 0.3 and 0.2 m and measures the Hausdorff distance to the study layout. | C and Dx are stable. D changes topology at h = 0.3 (8 lines instead of 12). This is reported as found, and the 8-line diagonal layout is studied as its own variant, Dd. |
| Mesh convergence with Richardson extrapolation | `verify/mesh_convergence.py` | Chooses the study mesh from the observed convergence rather than from a fixed rule. | h = 0.6 to 0.2 m for A, C, D. From h = 0.4 to 0.2 the change is at most 0.18 % (D) and 0.02 % or less for A and C; D converges unevenly (observed order 0.3), so its Richardson estimate is not relied on. Study mesh h = 0.4. |
| Rib effectiveness at equal volume | `verify/rib_effectiveness.py` | The ribbed variants buckle below the uniform shell, so this checks whether the ribs work at all. | The bare 50 mm shell carrying the V60 weight gives 221; C and D give 250 and 238, so the ribs add 7 to 13 %. Tripling the rib depth lowers lambda for D by the ratio of the added weight (0.757 against 0.75), and the mode-1 peak lies in the panels (peak on a rib node 0.37 to 0.41 of the shell peak). The governing mode is panel buckling between the ribs. |
| Beam-shell coupling on the curved shell | `verify/rib_coupling_curved.py` | The flat-slab check cannot show an error in rib orientation or offset sign on a curved surface. | Point load on a rib of D: deflection 0.468, 0.258 and 0.168 mm for rib depths 1, 209 and 800 mm. The ribs are coupled. |
| Mesh quality gate | `src/ribshell/mesh.py` | Rejects any mesh with min sin(corner angle) < 0.2. This keeps degenerate quads, which give a nonpositive Jacobian, out of the study. | Active for all runs. Earlier results are in `results/superseded/`. |
| Limit-point versus solver-stop classification | `study/s04_gnia.py` | CalculiX has no arc-length method, so each GNIA stop is classified by its end secant stiffness. A solver stop is never reported as a buckling load. | Applied to every GNIA row. |
| Both imperfection signs | `study/s04_gnia.py` | Runs the mode-1 imperfection in both signs. The weaker sign governs. | Applied to every GNIA row. |
| Nonlinear check at the study mesh | `study/s04_gnia.py`, `results/gnia_h0.4.txt` | h = 0.4 m, load control to 1.2 λ_LBA, amplitudes L/300, L/500 and L/1000, both signs. The recorded value is the weaker sign, and for these rows that sign is a limit point. | Finished 19 Nov 2026. Governing λ_GNIA: A 118.20 at L/500, C 91.30 at L/500, D 89.77 at L/300, E 135.96 at L/500. |
| Tool environment fix | `src/cecore/env.py`, `tests/test_env.py` | `setup_ccx` grew PATH on every call and crashed long runs (0xC0000017). It now adds only the missing entries, and a test checks this. | Fixed; the test passes. |
| Export round-trip check | `study/s07_exports.py`, `exports/ifc_check.txt` | Each IFC4 file is re-opened, validated against the schema, and its tessellated solid volume is compared with the model volume, so the BIM files carry the analysed geometry and not a drawing of it. | All six files: 0 schema errors; volume within 0.11 % of the model (A -0.02 %, B -0.11 %, C -0.08 %, D -0.02 %, Dd +0.03 %, Dx +0.10 %). The DXF passes the ezdxf audit. |
| One-command reproduction | `reproduce.py` | Runs the checks, the study, the exports and the figures in order, so a reader can rebuild every number from a fresh clone. | `--quick` rebuilds figures and exports from `results/` in a few minutes. |

## Choices to be confirmed by the author

| Choice | Where | Value and basis |
| --- | --- | --- |
| Stress-line tracing settings | `src/ribshell/layouts.py` | Anisotropy cut-off 0.02 for C. For D and Dx, lines stop below 15 % of the peak moment, where the principal directions are numerical noise. C and D are traced on the D4 group-averaged field with whole orbits. *Author to confirm.* |
| Edge clearance of ribs | `src/ribshell/layouts.py` | Lines whose median distance to the pinned boundary is below 0.5 m are dropped: the edge is already restrained, so such a rib only uses volume. *Author to confirm.* |
| Reference imperfection amplitude | `study/s04_gnia.py` | L/500 (L = 12 m plan span) as the reference, with L/300 and L/1000 as the bounds. To be checked against the code basis adopted for the report. *Author to confirm.* |
| Two-stage optimisation of variant E | `study/s03_optimise_E.py` | Stage 1: D4-symmetric, 6 variables, load directions 0° and 45°. Stage 2: free, 25 variables, 8 load directions, started from the stage-1 optimum. *Author to confirm.* |
| Study mesh size | `verify/mesh_convergence.py` | h = 0.4 m (element edge about 0.2 m after quad subdivision). h = 0.6 meets the scope's 2 % rule too, but D is 1.3 % stiff there against 0.2 % or less for A, which would favour the ribbed variants. *Author to confirm.* |
| Carbon factor | `study/s05_volume_sweep.py` | 0.120 kgCO2e/kg A1-A3 (band 0.063 to 0.149): IStructE (2022), *How to calculate embodied carbon*, 2nd ed., Table 2.3, UK in-situ C32/40 with 25 % GGBS, from the ICE database (2020). The table has no C30/37 row, so the next class up is used, which is on the high side. Reinforcement is excluded. *Author to confirm.* |
| Imperfect-model rib orientation | `src/ribshell/gnia.py`, `src/ribshell/fe.py` | The rib normals of the imperfect model come from the perfect surface. n1 is formed from the displaced tangent, so the section stays perpendicular to the rib; the tilt of n2 is about w0 / (half-wave) < 0.03 rad and is neglected. Stated as a limitation. |
| Vendoring `cecore` | `src/cecore/` | `cecore` (commit 8202eb8) is copied into this repository so that a fresh clone installs with `pip install -e .`; ShellBench is installed from its v1.0.0 DOI release. *Author to confirm.* |

## Record keeping

- The planning brief written before the study held private planning notes. It was replaced by the public scope document `docs/scope.md`, which keeps the research question, hypothesis, model, variants, methods, checks and outputs word for word, with only the timeline and the planning notes taken out. The brief was removed from the repository history so the public history matches the public scope.
