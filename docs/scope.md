# Study scope: stress-line ribs vs optimised thickness

**Working title:** How much buckling capacity does a buildable rib layout give up? Stress-line ribs vs optimised thickness in a form-found concrete shell under seismic sway, at equal embodied carbon

Author: Syed Nafiz Reja Hasib. Scope fixed 10 Nov 2026 before any results were computed; the hypothesis below is stated as it was then.

---

## 1. Design of the study

| Issue | How the study handles it |
|---|---|
| "Ribs beat uniform thickening at equal volume" is already known (bending stiffness grows with t³). Placing ribs along stress lines has also been published (Li et al. 2017, CGF; Tam & Mueller 2017). | The study now measures how close a buildable rib layout gets to an **optimised** material distribution (the upper bound), and whether the field you trace the ribs from matters. |
| Earthquake load is reversible and comes from any direction. Ribs traced for +x sway can be weak for y sway. | Every layout is scored by **λ_min over load directions** (0°, 45°, 90°, both signs). Ribs are traced from an envelope field of the x and y sway cases. |
| The volume comparisons didn't match (the 90 mm shell matched nothing). | Every variant is built at the same volume V, swept over V/V₆₀ = 1.0 to 1.5. Equal volume means equal A1–A3 carbon, so the main result is λ_cr vs tCO₂e. |
| Linear buckling (LBA) alone overstates the capacity of a dome. | Run a geometrically nonlinear check with the first-mode imperfection (GNIA) on three of the variants and report the knockdown. |

## 2. Research question

At equal concrete volume, how much of the lateral-load buckling capacity of an optimised variable-thickness shell does a rib layout traced from principal-moment lines recover? And does tracing from the seismic case instead of the gravity case change the answer?

**Hypothesis (falsifiable):** ribs traced from the seismic envelope (D) recover at least 70% of the gain that E achieves over A, while B and C recover clearly less. If B ends up close to D, the conclusion is that rib layout matters little, and that result is reported as it is.

## 3. Model

- **Geometry:** 12 × 12 m in plan, rise 2.4 m (span/5), pinned along all four edges. Form-found with your own force density code in NumPy (load: self-weight).
- **Material:** C30/37 concrete, E = 33 GPa, ν = 0.2, ρ = 2500 kg/m³. Linear elastic, uncracked (state this as a limitation).
- **Loads:** self-weight G plus an equivalent static lateral acceleration a_h = 0.20 g on the mass (sensitivity at 0.15 and 0.30 g). The lateral case is the one that drives buckling.
- **Elements:** S8R shells (CalculiX). Ribs are B32R beams with an offset, sharing nodes with the shell along curves embedded in the Gmsh surface (`gmsh.model.mesh.embed`). Cross-check one case with thickened shell strips, because beam–shell coupling in CalculiX has quirks.
- **Carbon:** A1–A3 only, one factor per m³ of concrete (cite the source; a reinforcement allowance is optional). Because every variant has equal volume, carbon is fixed within each set.

## 4. Variants (all at the same V)

| ID | Layout | What it tests |
|---|---|---|
| A | Uniform thickness t = V/area | Baseline |
| B | Orthogonal rib grid, about 1.5 m spacing | Ribs without mechanics-based placement |
| C | Ribs along gravity membrane-force lines | The classical "stress line" approach |
| D | Ribs along principal-moment lines of the seismic envelope | The proposed layout |
| E | Optimised thickness field (zones) | Upper bound, hard to build |

For B to D: base shell 50 mm, rib depth set so the total volume equals V. Rib width is fixed at 150 mm.

## 5. Methods

1. **Form-finding.** Force density method (NumPy, sparse solve). Check the residual force at each node and that the bending moment under gravity is close to zero.
2. **Linear static analysis.** Run G + a_h in the x and y directions. Export nodal membrane forces (n11, n22, n12) and moments (m11, m22, m12).
3. **Envelope field.** For each node, combine the moment tensors of the x and y cases into a symmetric envelope tensor (start with the sum of the absolute principal moments in tensor form; document the choice). Compare it with tracing from x only.
4. **Streamline tracer (own code).**
   - Interpolate the tensor field in plan view, since the shell is a height field z(x, y).
   - Integrate with RK4. At every step, flip the direction sign to match the previous step (tensor lines have no orientation).
   - Stop at boundaries and near degenerate points (where |σ1 − σ2| is small).
   - Seed with evenly spaced streamlines (Jobard–Lefer method), d_sep about 1.5 m.
   - Smooth and resample the curves, then embed them in the Gmsh surface.
   - Reference to compare against: Wang et al. 2026, *Topology-aware stress analysis in shell structures* (PSLshell code).
5. **Buckling (LBA).** CalculiX `*BUCKLE`, first 10 modes, λ on the combination G + a_h. Also report λ for gravity alone. Check how CalculiX handles a preload step in a buckle analysis before trusting the definition of λ.
6. **Optimisation (E).** 16–25 thickness zones, with mirror symmetry where the load envelope allows it. Objective: maximise λ_min over load directions. Constraints: equal volume and t between 40 and 150 mm. Use CMA-ES (`pycma`, free) with CalculiX in the loop; run it overnight. Watch for mode switching and repeated eigenvalues (Townsend & Kim 2019 explain why they break gradient methods).
7. **Nonlinear check (GNIA).** For A, D and E: scale the first buckling mode to an imperfection amplitude w₀ = L/500 (sensitivity at L/300 and L/1000), run a NLGEOM static analysis with increasing load until the limit point, and report knockdown = λ_GNIA / λ_LBA.
8. **Fabrication metric.** Total rib length, number of rib intersections and minimum rib curvature radius (for B to D). For E, the thickness range and gradient. This makes the "price of buildability" concrete.

## 6. Verification (do this first, before any study results)

| Check | Target | Pass criterion |
|---|---|---|
| Scordelis-Lo roof, S8R | 0.3024 ft vertical displacement at the free-edge midpoint | Within 2% after mesh convergence |
| Spherical cap buckling (classical/Zoelly) | p_cr = 2E t² / (R² √(3(1−ν²))) | LBA within 5% on a deep cap |
| Ribbed plate | Hand T-beam stiffness | Within 5% |
| Mesh convergence of λ_cr | 3 meshes | Change below 2% |
| Force density check | Bending under gravity | m ≈ 0, residual below 1e-6 |

## 7. Outputs

- Main figure: λ_min vs V (and tCO₂e) for A to E, plus a second axis or marker size for the fabrication metric.
- Buckling mode plots and a short animation of mode shapes (rendered with PyVista).
- Public repo: `form/`, `mesh/`, `fe/`, `trace/`, `opt/`, `verify/`, `figures/`, with one script that reproduces every figure.
- 8–10 page technical report (later a preprint). Journal options to decide once results are in: *Journal of the IASS*, *Structures*, *Engineering Structures*.


## 10. Room for improvement (only after the core is done)

- Corner-supported shell, which is much more sensitive to sway
- Response spectrum or time-history check instead of equivalent static load
- Cracked or nonlinear concrete
- Gradient-based thickness optimisation with adjoint sensitivities
- Optimising rib depth along each rib
- Second geometry (barrel vault)
- Test on a 3D-printed scale model

## 11. Out of scope (stated in the report)

Adjoint topology optimisation of the full rib network, reinforcement design, connections, construction-stage loads, and code compliance checks.

## Key references to read first

1. Li et al. (2017) Rib-reinforced shell structure. *Computer Graphics Forum* 36(7). doi:10.1111/cgf.13268
2. Tam & Mueller (2017) Additive manufacturing along principal stress lines. *3D Printing and Additive Manufacturing* 4(2). doi:10.1089/3dp.2017.0001
3. Oval, Rippmann, Mesnil, Van Mele, Baverel, Block (2019) Feature-based topology finding of patterns for shell structures. *Automation in Construction*. doi:10.1016/j.autcon.2019.02.008
4. Mesnil et al. (2025) Thickness optimisation in 3D printed concrete structures. *Automation in Construction*. doi:10.1016/j.autcon.2025.106076
5. Townsend & Kim (2019) A level set topology optimization method for the buckling of shell structures. *SMO* 60(5). doi:10.1007/s00158-019-02374-9
6. Wang, Liu, Wu, Westermann (2026) Topology-aware stress analysis in shell structures. Code: github.com/PSLer/PSLshell
