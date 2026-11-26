# Verified facts for the shell-buckling preprint

Source of truth for an 8–10 page computational note on the closed P2 study. Every factor below was read from a file in this repository at commit `2025f66` (20 Nov 2026, `main`). No CalculiX job was rerun for this sheet. Derived quantities are marked **derived**. Items the repository still marks for the author are marked **author to confirm**.

There is no `paper/REPO_MAP.txt`, no `docs/handoff.md`, and no `LOG.md`. The study record is `report/report.md`, `docs/scope.md` (fixed 10 Nov 2026, before the results), and `docs/decisions.md`.

Public repo: `syednafizreja/stress-line-ribs-shell-buckling`. Local path: `ComputationalEngineering/P2_ribs_vs_thickness`. Solver: CalculiX 2.23 (`verify/out/ccx_capabilities.txt`).

## 1. Geometry, loading, and equal-volume design space

Fixed in `src/ribshell/constants.py` and `docs/scope.md` section 3.

| Quantity | Value | Where |
|---|---|---|
| Plan | 12 m × 12 m square | `SPAN = 12.0` |
| Rise | 2.4 m (span/5), reached exactly | `RISE = 2.4`; `verify/out/fdm_residual.txt` prints `rise = 2.400000 m` |
| Form-finding | Force-density network, 60 × 60, self-weight, then a spline | `q_x = q_y = 6503.02 N/m`. Max relative vertical residual `5.583e-13` (limit `1e-6`). Max \|spline − network z\| `1.776e-15 m` |
| Supports | Translations fixed on every plan-boundary node. Rotations free | `*BOUNDARY` on set `EDGE`, degrees of freedom 1–3 (`src/ribshell/fe.py`) |
| Material | Linear-elastic, uncracked, called C30/37 | `E = 33 GPa`, `ν = 0.2`, `ρ = 2500 kg/m³` |
| Load | Self-weight plus a horizontal body force | `*DLOAD, GRAV`. Gravity is `g = 9.81 m/s²` in −z. Lateral acceleration is `a_h · g` with `a_h = 0.20`, direction `(cos θ, sin θ, 0)`. Sensitivity runs use `0.15` and `0.30`. No other live load |
| What λ scales | Self-weight and the lateral acceleration together. No preload | Deck kind `"combined"` in `fe.py`. `"lateral"` (gravity as a static preload, buckle step `*STEP, PERTURBATION`) exists in the code and is not the registered factor |

`λ_c` is the lowest positive buckling factor of that combination, over `θ = 0°, 45°, …, 315°` (eight directions, ten modes). `λ_g` is the same factor with self-weight alone. Variant E’s symmetric search scores only `0°` and `45°`, which is the whole orbit under D4 symmetry.

Registered volume: a uniform 60 mm shell,

`V60 = 9.669152495733323 m³`

in `results/variants_h0.4_r1.json`. **Derived** mid-surface area `V60 / 0.060 = 161.1525416 m²`.

Designs B–Dx use a 50 mm shell and ribs 150 mm wide. Depth is constant along each rib and is set by `src/ribshell/layouts.py`:

`d = (V − A_s t_b) / (w (L − k w))`

with `t_b = 50 mm`, `w = 150 mm`, `L` the rib length on the surface, and `k` the number of rib crossings. The `k w` term stops the crossing footprint being counted twice in the target volume. The finite-element volume sum in `Model.volume` does not subtract that footprint, so a layout with crossings comes out heavy (section 6).

The seven layouts at `V = V60`, mesh size `h = 0.4 m`:

| ID | Definition in the code | Ribs | Length | Crossings | Depth | Min. radius | Meshed volume |
|---|---|---:|---:|---:|---:|---:|---:|
| A | Uniform shell, 60 mm | 0 | 0 | 0 | — | — | 9.669 m³ |
| B | Orthogonal grid, spacing 1.5 m, both directions, lines at 1.5, 3.0, …, 10.5 m. No stress field | 14 | 179.469 m | 49 | 62.417 mm | straight | 9.738 m³ |
| C | Both families of gravity membrane-force lines. Anisotropy cut-off 0.02. D4 group average, whole orbits. Lines whose median distance to the pinned edge is below 0.5 m are dropped | 8 | 68.165 m | 0 | 157.604 mm | 1.818 m | 9.669 m³ |
| D | Both families of principal-moment lines of the seismic envelope `M_env = Σ_{±x,±y} Σ_i \|m_i\| v_i v_iᵀ`. Lines stop below 15% of the peak moment. Same D4 tracing and 0.5 m edge clearance | 12 | 51.443 m | 0 | 208.823 mm | 4.228 m | 9.669 m³ |
| Dd | The 8-line diagonal trace of the same moment field on the `h = 0.3 m` mesh, analysed at `h = 0.4 m`. Not in the registered hypothesis. D changes from 12 lines at `h = 0.4 m` to 8 lines at `h = 0.3 m` | 8 | 26.847 m | 0 | 400.176 mm | 4.519 m | 9.669 m³ |
| Dx | Moment lines of the +x case only (`abs_tensor` of `m`), mirror symmetry about `y = L/2` | 5 | 21.270 m | 0 | 505.083 mm | 3.544 m | 9.669 m³ |
| E | Shell thickness on a 5 × 5 zone grid, 40 mm to 150 mm, volume held at V60. Reported result is the D4 search: 6 variables, directions 0° and 45°, CMA-ES, 310 evaluations, stop `maxfevals = 300` | — | — | — | thickness 40.000–80.629 mm | — | 9.669 m³ |

Polylines: `results/layouts.json`. Depths and volumes: `results/variants_h0.4_r1.json`. E thicknesses: `results/E_best_d4_h0.4.json`, key `zone_thickness_m` (25 values, index `ix * 5 + iy`). A second stage, 25 free zones, eight directions, started from the D4 optimum, reaches objective `428.8755` and does not beat `429.0206`. It was run only at V60 (`results/E_best_free_h0.4.json`).

Tracing cut-offs, the 0.5 m edge clearance, the two-stage optimiser, and the mesh size are **author to confirm** (`docs/decisions.md`).

## 2. Numerical formulation and analytical verification

### Elements and mesh

- Shell: CalculiX `S8R`. Nodes of each quad are counter-clockwise in plan, so the normal points up. `src/ribshell/fields.py` states that CalculiX expands each S8R to a C3D20R through the thickness and integrates it with 2 × 2 × 2 Gauss points. Stress resultants used for tracing are recovered from those eight points.
- Ribs: `B32R`, rectangular section, one element set per beam because the section axis changes along a curve. Positive `OFFSET2 = (t_shell/2 + depth/2) / depth` hangs the rib below the mid-surface. Offset sign was checked in `verify/out/ccx_capabilities.txt` (T2).
- Study size `h = 0.4 m`. `docs/decisions.md` states that this is an element edge of about 0.2 m after quad subdivision. The saved mode mesh agrees: design A is 60 × 60 quadrilaterals (3600 elements, plan edge 0.2 m, 3721 corner nodes in `results/modes_h0.4.npz`). The npz stores corners only (`quads` has shape `(n, 4)`), not the eight S8R nodes. Embedded ribs remesh the shell: C has 4520 quads, D has 5512.
- Linear buckling deck: `*STEP` then `*BUCKLE` with 10 modes, then both gravity lines. No `SOLVER=` override. The factor used is the lowest positive eigenvalue. This is not a preload-and-perturbation buckle.
- Imperfect path: `*STATIC, NLGEOM`, load control, initial increment 0.05, minimum increment `1e-5`, at most 1000 increments, load ramped to `1.2 λ_LBA`. CalculiX 2.23 does not accept `*STATIC, RIKS`: `verify/out/ccx_capabilities.txt` (T3) records the parameter as unrecognised. A stop is a limit point only when the secant stiffness over the last 10% of the load is below 10% of the initial secant stiffness (`gnia.limit_kind`). A solver stop is not used as a buckling load.

### Benchmarks, copied from `verify/out/`

| Check | Reference | Finite element | Error | File |
|---|---|---|---|---|
| Scordelis–Lo roof, S8R, ShellBench model, same ccx build. Target is the free-edge midpoint vertical displacement 0.3024 (MacNeal & Harder 1985, as cited in `verify/scordelis_lo.py`) | 0.3024 | 8×8: 0.30325; 16×16: 0.30196; 24×24: 0.30192; 32×32: 0.30192 | finest mesh −0.16% (limit 2%). Pass | `verify/out/scordelis_lo.txt` |
| Spherical cap, Zoelly. `p_cr = 2 E (t/R)² / sqrt(3(1−ν²))`. Hemisphere `R = 5.0 m`, `t = 50 mm`, `E = 33 GPa`, `ν = 0.2`, pinned equator, uniform pressure. Mesh 27972 S8R, 83918 nodes, `h = 0.2 m` | 3.8891 MPa | `λ_1` gives 3.9224 MPa | +0.86% (script limit 5%). Pass | `verify/out/zoelly_cap.txt` |
| Ribbed slab, 6 m, 50 mm plate, eight ribs 150 × 250 mm | T-beam midspan 1.5242 mm | B32R ribs, inner four, 1.4883 mm, ratio 0.9764. Thickened S8R strips, ratio 0.9991. Beam/strip 0.9773 | within 5%. Pass | `verify/out/ribbed_strip.txt` |
| Mesh, `λ` at 0°, lowest positive of 10 modes, `h = 0.6, 0.4, 0.3, 0.2 m` | change from 0.4 m to 0.2 m | A `λ_c0`: 341.21 → 341.16 (−0.015%). D `λ_c0`: 235.11 → 234.69 (−0.179%). A `λ_g`: 369.77 → 369.72 | A and C change by 0.02% or less. D’s observed order on `λ_c0` is 0.32, so the Richardson estimate for D is not used. Study mesh stays `h = 0.4 m`. **Author to confirm** | `verify/out/mesh_convergence.txt` |
| Curved rib coupling, 10 kN on rib 9 of D | deflection must fall as the rib deepens | `w = 0.4676, 0.3989, 0.3319, 0.2578, 0.2069, 0.1684 mm` at depths 1, 50, 100, 209, 400, 800 mm | ribs are coupled | `verify/out/rib_coupling_curved.txt` |

The scope also asked for a near-zero bending moment under gravity on the form-found shell. No such moment is stored. The force-density file stores only the vertical residual and the spline error.

## 3. Buckling and imperfection-sensitivity results

### Linear factors at V60, `h = 0.4 m`, `a_h = 0.20 g`

Unrounded `lambda_c` from `results/variants_h0.4_r1.json`. E from `results/E_best_d4_h0.4.json`, key `objective`. Percentages are **derived** from those values against A’s `341.2093`. The landing-page rounding of the same ratios is the last column.

| Design | `λ_g` | `λ_c` | Worst θ | Against A, unrounded | Landing-page rounding |
|---|---:|---:|---:|---:|---:|
| A | 369.7709 | 341.2093 | 180° | baseline | baseline |
| B | 428.3312 | 395.9731 | 0° | +16.0499% | +16.05% |
| C | 268.5123 | 249.1289 | 180° | −26.9865% | −26.99% |
| D | 259.9983 | 235.0948 | 180° | −31.0995% | −31.10% |
| Dd | 290.3324 | 268.9580 | 180° | −21.1751% | −21.17% |
| Dx | 211.7727 | 193.1005 | 0° | −43.4070% | −43.41% |
| E | not stored | 429.0206 | min of 0° and 45° is 0°, both 429.0206 at 0° and 438.0098 at 45° | +25.7353% | +25.74% |

`E − A = 87.8113`. **Derived** recovery of that gain: B recovers `0.6237` (62.37%). D is below A, so it recovers none (`R_D = −1.208` in `results/sweep_h0.4.json` at `V/V60 = 1`). C is also below A. The registered hypothesis (D recovers at least 70% of E’s gain, and B and C recover less than D) is not supported.

A, B, C, D, and Dd stay within about 0.01% across a square’s symmetry orbit. Dx, traced from +x only, spreads from 193.1005 at 0° to 233.795 at 180°, which is 21.07% of `λ_c` inside one orbit (`results/variants_h0.4_r1.txt`).

Printed two-decimal table used on the landing page: A 341.21, B 395.97, C 249.13, D 235.09, Dd 268.96, Dx 193.10, E 429.0206.

### Mode-1 imperfect factors

Run for A, C, D, and E only (`study/s04_gnia.py`, `results/gnia_h0.4.json`). B, Dd, and Dx have no geometrically nonlinear result.

The first mode is scaled so its largest translation equals `w0`. The sign called `peak-down` moves that node downward. Both signs are run. `L = 12 m`, so the three amplitudes are `L/300 = 0.040 m`, `L/500 = 0.024 m`, and `L/1000 = 0.012 m`. The driver picks the lower of the 0° and 45° eigenvalues as the imperfection shape. For these four designs that angle is 0°. The linear factor stored beside D in the nonlinear file is therefore the 0° value `235.1121`, not the eight-direction minimum `235.0948`.

Governing value of each row: the weaker sign. For every governing row below, that sign is `peak-down` and the classifier says `limit point`. The other sign is often a solver stop (`figures/table_gnia.csv`).

| Design | `w0` | `λ_LBA` (0°) | `λ_GNIA` | `α = λ_GNIA / λ_LBA` |
|---|---|---:|---:|---:|
| A | L/300 | 341.209400 | 121.768722 | 0.356874 |
| A | L/500 | 341.209400 | 118.202975 | 0.346424 |
| A | L/1000 | 341.209400 | 148.136361 | 0.434151 |
| C | L/300 | 249.132400 | 105.569854 | 0.423750 |
| C | L/500 | 249.132400 | 91.295285 | 0.366453 |
| C | L/1000 | 249.132400 | 125.088671 | 0.502097 |
| D | L/300 | 235.112100 | 89.771960 | 0.381826 |
| D | L/500 | 235.112100 | 98.291689 | 0.418063 |
| D | L/1000 | 235.112100 | 118.371146 | 0.503467 |
| E | L/300 | 429.020600 | 140.405314 | 0.327269 |
| E | L/500 | 429.020600 | 135.958464 | 0.316904 |
| E | L/1000 | 429.020600 | 165.937120 | 0.386781 |

Minimum over the three amplitudes, weaker sign: A `118.202975` at L/500, C `91.295285` at L/500, D `89.771960` at L/300, E `135.958464` at L/500.

The published nonlinear comparison of E with A uses L/500 for both, which is also each design’s minimum. From the JSON,

`(135.958464 − 118.202975) / 118.202975 = +15.021%`.

The report’s `+15.03%` is the same ratio after rounding each factor to two decimals: `(135.96 − 118.20) / 118.20 = +15.025%`, which rounds to `+15.03%`. Quote one of these and say which. Do not quote `+15.03%` as if it were the unrounded ratio.

At that same L/500 amplitude, C is `(91.295285 − 118.202975) / 118.202975 = −22.76%` against A, and D is `(98.291689 − 118.202975) / 118.202975 = −16.84%` against A. D’s lowest nonlinear factor is the L/300 row, `89.771960`, which is `−24.05%` against A’s L/500 factor. The order E > A > C > D is unchanged. B, Dd, and Dx are unknown on this path.

`L/500` as the reference amplitude is **author to confirm** against the code basis used for design.

### Volume sweep and lateral acceleration (context, not the headline table)

`results/sweep_h0.4.txt`, `a_h = 0.20 g`, E re-optimised at each volume with the six-variable D4 search. C, D, Dd, and Dx stay below A at every ratio from 1.0 to 1.5, and their `λ_c` falls as the volume rises. B stays above A. From 1.1 to 1.3 times V60, B is above the recorded D4 optimum. The 25-zone search was not repeated off V60, so B above E at `V/V60 = 1.2` (`684.7` against `546.9`) is a comparison with the symmetric search only.

`results/ah_h0.4.txt`: the same V60 designs at 0.15 g, 0.20 g, and 0.30 g. E is the recorded D4 field, not re-optimised. The order is E, B, A, Dd, C, D, Dx at all three accelerations.

## 4. Mechanical mechanism

What the files support, and what they do not.

Designs B–Dx buy their ribs by thinning the shell from 60 mm to 50 mm. **Derived** from those two thicknesses: plate bending stiffness of an unreinforced panel scales as `t³`, so a 50 mm panel has `(50/60)³ = 0.5787` of the bending stiffness of A. Membrane stiffness scales as `t`, so the same panel has `5/6 = 0.8333` of A’s membrane stiffness. No modal membrane-energy or bending-energy split is stored, so a membrane-to-bending shift is not a measured result.

The stress-line layouts do not subdivide the bay. C, D, Dd, and Dx have zero crossings (`results/variants_h0.4_r1.json`). Their spare volume, about one sixth of V60, sits in short, deep ribs: D is 12 lines, 51.443 m, 208.8 mm deep. B is the opposite geometry: 14 lines, 179.469 m, 49 crossings, 1.5 m spacing, only 62.4 mm deep, and it is the rib layout that beats A (`+16.05%`).

A separate check at `h = 0.6 m` and `θ = 0°` only (`verify/rib_effectiveness.py`, `verify/out/rib_effectiveness.txt`) asks whether the ribs are working:

| Model | `λ` at 0° |
|---|---:|
| A, uniform 60 mm | 341.41 |
| Bare 50 mm shell, density scaled by 60/50 so the weight equals V60 | 221.23 |
| C at the equal-volume depth, 158 mm | 249.75 (1.129 × the bare shell) |
| D at the equal-volume depth, 209 mm | 237.68 (1.074 × the bare shell) |
| C with ribs three times deeper, volume not held, 473 mm | 227.56 (0.911 × that C) |
| D with ribs three times deeper, 626 mm | 179.92 (0.757 × that D) |

On the mode of that check, the largest displacement on a rib node is 0.370 (C) and 0.409 (D) of the shell peak, so the peak is in the shell. The mean rib-node amplitude is 0.992 (C) and 1.345 (D) times the mean shell amplitude, so D’s ribs are not nodal lines. Tripling D’s depth multiplies `λ` by 0.757. If the rib stiffness did not change and only the mass changed, the equal-volume shell (60 mm worth of material) against a shell with the rib volume tripled (80 mm worth) would give a factor `0.060/0.080 = 0.75`. The recorded 0.757 sits on that mass ratio. Extra rib depth did not buy buckling stiffness. `docs/decisions.md` calls the governing mode panel buckling between the ribs. That sentence matches the peak ratios and the mass-proportional drop. It does not match a picture in which the ribs are quiet nodal lines.

This check is not the study mesh and not the minimum over eight directions. It was not run for B, Dd, Dx, or E.

E keeps one continuous shell. At V60 the D4 optimum runs from 40.000 mm, the lower bound, to 80.629 mm (`zone_thickness_m`). Nine of the 25 zones sit on 40 mm. Zones are 2.4 m square. The largest step between zones that share an edge is 40.6 mm, a plan gradient of 0.0169 (`report/report.md` section 8, taken from the saved thicknesses). The search objective is `λ_c`. It raises the factor from the uniform shell’s 341.2094 (the `uniform_or_start_objective` stored in the same JSON) to 429.0206. That is a redistribution of the same volume through the shell thickness. It is not evidence about reinforcement, cracking, or a measured energy split.

## 5. Figure data mapping

Four figures are enough for the note. Plot from these files. Do not read `.frd` or `.dat` dumps. `study/s06_figures.py` already draws six figures from the same files; the mapping below is the preprint set.

### Figure 1 — seven designs at one volume

Bar chart of `λ_c`, with `λ_g` as a second series for A–Dx.

- `results/variants_h0.4_r1.json` → `variants.<ID>.lambda_c`, `lambda_g`, `worst_angle`, `rib_length_m`, `n_ribs`, `crossings`, `rib_depth_mm`, `volume_check_m3`.
- E is not in that file. Use `results/E_best_d4_h0.4.json` → `objective` (429.0206). There is no `λ_g` for E.
- Rounded copy, A–Dx only: `figures/table_variants_V60.csv` columns `variant, lambda_g, lambda_c, worst_angle, rib_length_m, n_ribs, crossings, rib_depth_mm`.
- Rib lines, if the bars are paired with a plan: `results/layouts.json`, keys `B`, `C`, `D`, `Dd`, `Dx`, each a list of plan polylines.

### Figure 2 — `λ_c` against volume

- `results/sweep_h0.4.json` → `rows[].vratio`, `rows[].V_m3`, `rows[].tCO2e.default` (also `low`, `high`), `rows[].lambda_c.{A,B,C,D,Dd,Dx,E}`, `rows[].recovery.{B,C,D,Dd,Dx}`, `rows[].fabrication.<rib ID>.{rib_length_m, crossings, rib_depth_mm, min_radius_m, n_ribs}`, `rows[].E_free` (present at V60 only, value 428.8755).
- Carbon factor `ecf_kgco2e_per_kg`: default 0.12, low 0.063, high 0.149 kgCO₂e/kg. **Author to confirm.** The source string in the file is IStructE 2022, *How to calculate embodied carbon*, 2nd ed., Table 2.3, UK in-situ C32/40 with 25% GGBS, from ICE 2020. The table has no C30/37 row. Reinforcement is not included. Top-axis tonnes are `V_m3 * 2500 * 0.12 / 1000`.
- Rounded copy: `figures/table_sweep.csv`.

### Figure 3 — imperfect paths

- `results/gnia_h0.4.json` → `rows[]` with `variant`, `w0`, `w0_m`, `sign`, `angle_deg`, `lambda_LBA`, `lambda_GNIA`, `alpha`, `kind`, and `history` as pairs `(λ, |u| at the mode-1 peak node in metres)`.
- Left panel: rows with `w0 == "L/500"`, both signs. Plot `|u|` in millimetres against `λ / lambda_LBA`. Draw solver-stop curves differently from limit points (`kind`).
- Right panel: for each of A, C, D, E and each of L/300, L/500, L/1000, the minimum `alpha` over the two signs.
- Rounded copy: `figures/table_gnia.csv` columns `variant, w0, sign, angle_deg, lambda_LBA, lambda_GNIA, alpha, kind`. The history is only in the JSON.
- Do not put B, Dd, or Dx on this figure.

### Figure 4 — where the mode moves, and where E puts the thickness

- Mode: `results/modes_h0.4.npz`. For `k` in `{A, C, D}`: `k_xyz` corner coordinates `(n, 3)`, `k_quads` corner connectivity `(m, 4)`, `k_amp` nodal `|u| / max`, `k_factor`, `k_angle`. Stored factors are A 341.2093 at 180°, C 249.1289 at 180°, D 235.0948 at 180°. These are the eight-direction minima, not the 0° shapes used for the imperfect path. Overlay polylines from `results/layouts.json` for C and D. E, B, Dd, and Dx are not in this npz.
- Thickness: `results/E_best_d4_h0.4.json` → `zone_thickness_m`, reshape `(5, 5)` with index `ix * 5 + iy`, times 1000 for millimetres. Optional second grid from `results/E_best_free_h0.4.json`, same key. `s06_figures.py` transposes for `imshow` so that rows run with `y`.

## 6. Candid limitations in the code

1. The material is linear-elastic and uncracked. Cracking, reinforcement, creep, shrinkage, and long-term modulus are absent (`docs/scope.md` sections 3 and 11; `report/report.md` section 11). The carbon row is concrete A1–A3 only, and the factor is the C32/40 row rather than C30/37. **Author to confirm.**
2. The lateral load is a constant acceleration of 0.2 g on the mass, scaled with self-weight. It is not a response spectrum and not a ground-motion history. On the uniform shell the sway cuts `λ` by about 8% at 0.2 g (`λ_c / λ_g = 0.923`), so this load case barely drives buckling. See the discrepancy below.
3. The imperfect analysis is one shape, mode 1, at one direction (the lower of 0° and 45°), at three amplitudes and two signs, and only for A, C, D, and E. Load control cannot pass a limit point. There is no Riks path. B, Dd, and Dx have no imperfect factor. Rib normals on the imperfect mesh are taken from the perfect surface (`docs/decisions.md`).
4. One geometry, constant rib depth, no adjoint layout search, no construction stages, no code check, and no physical model. The corner-supported shell was screened and kept out of the comparison because its sway ratio stayed above 0.91 (`report/report.md` section 9). The scope’s “bending under gravity is close to zero” check has no stored number.

## DISCREPANCIES TO RESOLVE

1. **`+15.03%` is a rounded-table ratio.** Unrounded L/500 factors give `(135.958464 − 118.202975) / 118.202975 = +15.021%`. The report and the README get `+15.03%` from `(135.96 − 118.20) / 118.20`. Use `+15.02%` from the JSON, or show the two-decimal arithmetic. Do not treat `15.03` as the file ratio.
2. **D has two linear factors, 0.017 apart.** Eight-direction `λ_c` is 235.0948 (180°). The nonlinear file and the 0° mesh-convergence column use 235.1121. The landing page prints 235.09. Both are real. They must not be mixed inside one percentage.
3. **B is not at V60 in the mesh.** Target volume 9.669152 m³. `volume_check_m3` for B is 9.737968 m³, +0.712%. The depth formula subtracts the 49 crossings; `Model.volume` adds `width × depth × length` with no crossing correction. A, C, D, Dd, Dx, and E match V60 to the stored digits. The `+16.05%` for B is slightly on the heavy side of the equal-volume claim.
4. **The scope says the lateral case drives buckling. The probe does not.** `docs/scope.md` section 3. `verify/out/sway_probe.txt` and `report/report.md` section 9: for the edge-pinned 60 mm shell, `λ_c / λ_g` is 0.943, 0.923, 0.884, and 0.849 at 0.15 g, 0.20 g, 0.30 g, and 0.40 g. The report already states this. The preprint should not repeat the scope sentence.
5. **Imperfect factors do not exist for all seven layouts.** The scope asked for A, D, and E. The run added C. B, Dd, and Dx were not run. D’s governing nonlinear row is L/300 (`89.77`), while A, C, and E govern at L/500. The README hero pairs 118.20, 91.30, 89.77, and 135.96. That is the minimum over amplitude, not one amplitude.
6. **No Riks, and the registered buckle is not a perturbation step.** `*STATIC, RIKS` is rejected by this CalculiX 2.23. `λ_c` comes from `*BUCKLE` on `G + 0.2 g` with no preload. The perturbation deck (`kind = "lateral"`) is a different definition and is not the table in section 3.
7. **Gravity bending moment was required and not filed.** Scope sections 5 and 6 ask for a near-zero moment under self-weight. `verify/out/fdm_residual.txt` has the network residual only. `verify/out/check_fields.txt` is a flat-plate patch test, not this shell.
8. **`paper/REPO_MAP.txt` is not in the repository.** This sheet was taken from `report/report.md`, `docs/scope.md`, `docs/decisions.md`, `src/ribshell/`, `results/`, and `verify/out/`.
