"""Fast checks of the model building blocks. None of these needs CalculiX, so they run in CI."""
from __future__ import annotations

import numpy as np
import pytest

from cecore.tracer import symmetry_ops
from ribshell.constants import SPAN, T_MAX, T_MIN, T_REF
from ribshell.form import shell_surface
from ribshell.gnia import limit_kind
from ribshell.layouts import Layout, crossings, grid_layout, ribbed_model
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign, d4_orbits


@pytest.fixture(scope="module")
def surface():
    return shell_surface()


@pytest.fixture(scope="module")
def mesh(surface):
    return build_mesh(surface, 0.6)


def test_form_is_pinned_and_rises(surface):
    # pinned plan boundary at z = 0, crown at the scope's rise (span / 5)
    assert abs(float(surface.z(0.0, SPAN / 2))) < 1e-6
    assert float(surface.z(SPAN / 2, SPAN / 2)) == pytest.approx(SPAN / 5, rel=1e-3)


def test_mesh_passes_quality_gate(mesh):
    assert mesh.corner_quality().min() >= 0.2
    # curved area exceeds the plan area but not by much for a rise of span / 5
    assert SPAN**2 < mesh.area() < 1.3 * SPAN**2


def test_d4_orbits_of_5x5_zones():
    orb = d4_orbits(5)
    assert orb.max() + 1 == 6
    g = orb.reshape(5, 5)
    assert np.array_equal(g, g.T) and np.array_equal(g, g[::-1]) and np.array_equal(g, g[:, ::-1])


def test_symmetry_group_has_eight_elements():
    ops = symmetry_ops("D4")
    assert len(ops) == 8
    for a in ops:
        assert np.allclose(a @ a.T, np.eye(2))


@pytest.mark.parametrize("sym", [None, "d4"])
def test_zonal_design_keeps_volume_and_bounds(mesh, sym):
    v = mesh.area() * T_REF
    des = ZonalDesign(mesh, 5, v, symmetry=sym)
    rng = np.random.default_rng(1)
    for _ in range(5):
        tz = des.zone_thickness(rng.uniform(0.1, 3.0, des.n_vars))
        assert des.zone_area @ tz == pytest.approx(v, rel=1e-10)
        assert tz.min() >= T_MIN - 1e-12 and tz.max() <= T_MAX + 1e-12


def test_zonal_design_rejects_infeasible_volume(mesh):
    with pytest.raises(ValueError):
        ZonalDesign(mesh, 5, mesh.area() * (T_MIN / 2))


def test_grid_crossings():
    lay = grid_layout()
    nx = sum(1 for ln in lay.lines if abs(ln[0, 1] - ln[-1, 1]) < 1e-9)
    assert crossings(lay.lines) == nx * (len(lay.lines) - nx)


def test_ribbed_model_hits_target_volume(surface, mesh):
    v = mesh.area() * T_REF
    lay = Layout("X", [np.array([[1.0, 6.0], [11.0, 6.0]]), np.array([[6.0, 1.0], [6.0, 11.0]])])
    m, met = ribbed_model(lay, surface, v, h=0.6)
    assert met["crossings"] == 1
    assert met["volume_m3"] == pytest.approx(v, rel=1e-9)
    # the FE volume double-counts the crossing; the difference is exactly one rib-width square
    assert m.volume() - v == pytest.approx(met["rib_depth_mm"] / 1e3 * 0.15**2, rel=0.05)


def test_limit_kind():
    lam = np.linspace(0.1, 1.0, 10)
    assert limit_kind(list(zip(lam, lam))) == "solver stop"            # still linear
    u = np.concatenate([lam[:8], [1.5, 4.0]])
    lam2 = np.concatenate([lam[:8], [0.85, 0.86]])
    assert limit_kind(list(zip(lam2, u))) == "limit point"             # flat end
    assert limit_kind([(0.1, 0.1), (0.2, 0.2)]) == "undetermined"
