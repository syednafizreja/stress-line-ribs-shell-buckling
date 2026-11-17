# Project 2: stress-line ribs vs optimised thickness

Brief: `brief.md`. Status: `../STATUS.md`. Local only: no GitHub remote until the author decides to release.

| Folder | Holds |
|---|---|
| `form/` | form-finding wrapper around `shellbench.formfind` (12 x 12 m, rise 2.4 m, pinned edges) |
| `mesh/` | Gmsh surface + embedded rib curves (`gmsh.model.mesh.embed`) |
| `fe/` | CalculiX decks: static G + a_h, `*BUCKLE`, NLGEOM |
| `trace/` | P2-specific seeding; the tracer itself lives in `../cecore` (shared with P1) |
| `opt/` | variant E, CMA-ES (`cma`) over 16-25 thickness zones |
| `verify/` | checks that must pass first (brief section 6) and `ccx_capabilities.py` |
| `figures/` | one script per figure |

Environment: conda env `research` + `pip install -e ../cecore`. Run every script with
`~/miniconda3/envs/research/python.exe`. CalculiX is found through `cecore.env.setup_ccx()`.

Smoke tests: `python verify/smoke_env.py` and `python verify/ccx_capabilities.py`.
