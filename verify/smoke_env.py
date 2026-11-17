"""Environment smoke test for P1 and P2: every import works, ccx runs, gmsh meshes,
and ShellBench's force density form-finding gives the P2 geometry.

This is a set-up check, not the scope's force density verification (that one also
checks the residual and m ~ 0 under gravity on the FE model).
"""
import importlib

import numpy as np

for m in ("numpy", "scipy", "matplotlib", "gmsh", "meshio", "cma", "pymoo", "pyvista",
          "ifcopenshell", "ezdxf", "imageio",
          "shellbench.formfind", "shellbench.calculix", "shellbench.diaphragm",
          "shellbench.shellfe", "shellbench.generator", "shellbench.rsa", "shellbench.carbon",
          "cecore.env"):
    mod = importlib.import_module(m)
    print(f"ok  {m:24s} {getattr(mod, '__version__', '')}")

from cecore.env import ccx_version
print("ok  ccx", ccx_version())

import gmsh
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.addRectangle(0, 0, 0, 12, 12); gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
print("ok  gmsh nodes:", len(gmsh.model.mesh.getNodes()[0])); gmsh.finalize()

from shellbench.formfind import form_find
v = form_find(12.0, 12.0, 2.4, w=2500 * 9.81 * 0.06, nx=40, ny=40)
print(f"ok  form_find 12x12, rise {v.z.max():.3f} m, n_x = {v.q_x * v.dx / v.dy / 1e3:.1f} kN/m (plan-uniform load)")
assert np.isclose(v.z.max(), 2.4)
