"""Model exports and renders of the V = V60 designs at the study mesh.

1. Mode 1 under the governing sway direction for A, C, D and E (CalculiX *BUCKLE), saved to
   results/modes_h<h>.npz for fig5 and the renders.
2. IFC4 (exports/P2_<variant>.ifc): the shell as an IfcSlab (ROOF) with a closed triangulated
   solid (mid-surface offset by +/- t/2 along the surface normal; E uses the nodal mean of the
   zone thickness), each rib as an IfcMember box chain hanging below the shell exactly as in the
   FE model (top face on the shell soffit), material C30/37 and a property set with thickness,
   volume, A1-A3 carbon and the buckling factors. The script re-opens each file, checks the
   schema and compares the closed-solid volume with the FE volume (exports/ifc_check.txt).
3. DXF (exports/P2_layouts.dxf, metres): one layer per rib layout as 3D polylines on the shell
   surface (plan view in CAD = the layout), shell contours every 0.3 m, and the E zone grid with
   the optimised thickness of each zone.
4. PyVista renders (figures/render_*.png) and a mode-1 animation of D (figures/mode_D.gif).

Usage: python study/s07_exports.py [--h 0.4]
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np

from ribshell.constants import RHO, T_BASE_RIBBED, T_REF
from ribshell.fe import Model
from ribshell.form import shell_surface
from ribshell.gnia import lba_mode
from ribshell.layouts import Layout, ribbed_model
from ribshell.mesh import build_mesh
from ribshell.zonal import ZonalDesign

ROOT = Path(__file__).resolve().parents[1]
RES, FIG, EXP = ROOT / "results", ROOT / "figures", ROOT / "exports"
ECF = 0.120                       # kgCO2e/kg A1-A3, IStructE 2022 Table 2.3 (see s05)
RIBBED = ("B", "C", "D", "Dd", "Dx")


# ---------------------------------------------------------------------------- models
from matplotlib.colors import LinearSegmentedColormap
#: grey (no displacement) to deep blue: the undisplaced shell stays visible on the light page
MODE_CMAP = LinearSegmentedColormap.from_list("grey_blue", ["#d9d8d3", "#9cc3ea", "#2a78d6", "#0b2f6b"])

def build_models(h: float) -> tuple[object, dict[str, Model], dict[str, dict]]:
    s = shell_surface()
    mesh = build_mesh(s, h)
    v60 = mesh.area() * T_REF
    var = json.load(open(RES / f"variants_h{h:g}_r1.json"))["variants"]
    lay = json.load(open(RES / "layouts.json"))
    ms = {"A": Model(mesh, T_REF, surface=s)}
    meta = {"A": {**var["A"], "thickness": "60 mm uniform"}}
    for n in RIBBED:
        if n in lay and n in var:
            ms[n], met = ribbed_model(Layout(n, [np.array(x) for x in lay[n]]), s, v60, h=h)
            meta[n] = {**var[n], **met, "thickness": f"{T_BASE_RIBBED * 1e3:.0f} mm shell + ribs"}
    eb = RES / f"E_best_d4_h{h:g}.json"
    if eb.exists():
        e = json.load(open(eb))
        des = ZonalDesign(mesh, e["nz"], v60)
        tz = des.zone_thickness(np.asarray(e["zone_thickness_m"]))
        ms["E"] = Model(mesh, des.element_thickness(tz), surface=s)
        lc = {int(k): v for k, v in e["lambda_c_by_angle"].items()}
        meta["E"] = {"lambda_c": e["objective"], "worst_angle": min(lc, key=lc.get), "zone_thickness_m": list(tz),
                     "thickness": f"zonal {tz.min() * 1e3:.0f}-{tz.max() * 1e3:.0f} mm", "lambda_g": None}
    return s, ms, meta


def corner_mesh(m: Model) -> tuple[np.ndarray, np.ndarray, list[int], dict[int, int]]:
    """Quad corners only (the S8R mid-side nodes are not needed for a surface)."""
    used = sorted({n for c in m.mesh.quads.values() for n in c[:4]})
    idx = {n: i for i, n in enumerate(used)}
    xyz = np.array([m.mesh.nodes[n] for n in used])
    quads = np.array([[idx[n] for n in c[:4]] for c in m.mesh.quads.values()])
    return xyz, quads, used, idx


def node_thickness(m: Model, used: list[int]) -> np.ndarray:
    acc: dict[int, list[float]] = {}
    for e, c in m.mesh.quads.items():
        for n in c[:4]:
            acc.setdefault(n, []).append(m.t_of(e))
    return np.array([np.mean(acc[n]) for n in used])


def shell_solid(m: Model, s) -> tuple[np.ndarray, np.ndarray]:
    """Closed triangulated solid of the shell: top, bottom and the side faces on the boundary."""
    xyz, quads, used, _ = corner_mesh(m)
    nrm = s.normal(xyz[:, 0], xyz[:, 1])
    t = node_thickness(m, used)[:, None]
    top, bot = xyz + nrm * t / 2, xyz - nrm * t / 2
    nv = len(xyz)
    tris = []
    for q in quads:                                   # quads are counter-clockwise in plan
        tris += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        tris += [[q[0] + nv, q[2] + nv, q[1] + nv], [q[0] + nv, q[3] + nv, q[2] + nv]]
    edges: dict[tuple[int, int], int] = {}
    for q in quads:
        for a, b in zip(q, np.roll(q, -1)):
            key = (min(a, b), max(a, b))
            edges[key] = edges.get(key, 0) + 1
    for q in quads:
        for a, b in zip(q, np.roll(q, -1)):
            if edges[(min(a, b), max(a, b))] == 1:    # boundary edge, oriented a -> b
                tris += [[a, b + nv, b], [a, a + nv, b + nv]]
    pts = np.vstack([top, bot])
    tri = np.array(tris)
    if signed_volume(pts, tri) < 0:
        tri = tri[:, ::-1]
    return pts, tri


def rib_boxes(m: Model, s) -> list[tuple[np.ndarray, np.ndarray]]:
    """One closed box chain per rib: rib top face on the shell soffit, depth below (as in fe.py)."""
    shell_at = {n: e for e, c in m.mesh.quads.items() for n in c}
    out = []
    for r in m.ribs:
        P, T = [], []
        for a, mid, b in m.mesh.ribs[r.rib_id]:
            for p0, p1 in ((a, mid), (mid, b)):
                q = [np.array(m.mesh.nodes[p0]), np.array(m.mesh.nodes[p1])]
                tv = (q[1] - q[0]) / np.linalg.norm(q[1] - q[0])
                corners = []
                for pt, nid in zip(q, (p0, p1)):
                    n2 = s.normal(pt[0], pt[1])
                    n1 = np.cross(n2, tv)
                    n1 /= np.linalg.norm(n1)
                    top = pt - n2 * m.t_of(shell_at[nid]) / 2
                    for dz in (0.0, r.depth):
                        for sy in (-1, 1):
                            corners.append(top - n2 * dz + n1 * sy * r.width / 2)
                o = len(P)
                P += corners          # 0..3 at p0: (top,-), (top,+), (bot,-), (bot,+); 4..7 at p1
                f = [[0, 1, 3], [0, 3, 2], [4, 7, 5], [4, 6, 7], [0, 4, 5], [0, 5, 1],
                     [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]]
                box = np.array(f) + o
                pts = np.array(P)
                if signed_volume(pts, box) < 0:
                    box = box[:, ::-1]
                T += box.tolist()
        out.append((np.array(P), np.array(T)))
    return out


def signed_volume(pts: np.ndarray, tri: np.ndarray) -> float:
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)


# ---------------------------------------------------------------------------- modes
def compute_modes(ms: dict[str, Model], meta: dict, h: float) -> dict:
    out = {}
    with tempfile.TemporaryDirectory() as td:
        for k in ("A", "C", "D", "E"):
            if k not in ms:
                continue
            ang = int(meta[k]["worst_angle"])
            r = lba_mode(ms[k], Path(td), f"mode_{k}", ang, threads=8)
            xyz, quads, used, _ = corner_mesh(ms[k])
            u = np.array([r.mode.get(n, np.zeros(3)) for n in used])
            amp = np.linalg.norm(u, axis=1)
            out[k] = {"xyz": xyz, "quads": quads, "u": u / amp.max(), "amp": amp / amp.max(),
                      "factor": r.factor, "angle": ang}
            print(f"mode {k}: lambda {r.factor:.2f} at {ang} deg", flush=True)
    np.savez_compressed(RES / f"modes_h{h:g}.npz", **{f"{k}_{f}": v[f] for k, v in out.items() for f in v})
    return out


# ---------------------------------------------------------------------------- IFC
def write_ifc(name: str, m: Model, s, meta: dict, h: float) -> dict:
    import ifcopenshell
    import ifcopenshell.api as api
    import ifcopenshell.api.aggregate
    import ifcopenshell.api.context
    import ifcopenshell.api.geometry
    import ifcopenshell.api.material
    import ifcopenshell.api.project
    import ifcopenshell.api.pset
    import ifcopenshell.api.root
    import ifcopenshell.api.spatial

    f = ifcopenshell.api.project.create_file(version="IFC4")
    proj = api.root.create_entity(f, ifc_class="IfcProject", name="Stress-line ribs vs optimised thickness")
    units = [f.createIfcSIUnit(None, "LENGTHUNIT", None, "METRE"),
             f.createIfcSIUnit(None, "AREAUNIT", None, "SQUARE_METRE"),
             f.createIfcSIUnit(None, "VOLUMEUNIT", None, "CUBIC_METRE"),
             f.createIfcSIUnit(None, "MASSUNIT", "KILO", "GRAM"),
             f.createIfcSIUnit(None, "PLANEANGLEUNIT", None, "RADIAN")]
    proj.UnitsInContext = f.createIfcUnitAssignment(units)
    model = api.context.add_context(f, context_type="Model")
    body = api.context.add_context(f, context_type="Model", context_identifier="Body",
                                   target_view="MODEL_VIEW", parent=model)
    site = api.root.create_entity(f, ifc_class="IfcSite", name="Study site")
    bldg = api.root.create_entity(f, ifc_class="IfcBuilding", name=f"Variant {name}")
    storey = api.root.create_entity(f, ifc_class="IfcBuildingStorey", name="Roof")
    api.aggregate.assign_object(f, products=[site], relating_object=proj)
    api.aggregate.assign_object(f, products=[bldg], relating_object=site)
    api.aggregate.assign_object(f, products=[storey], relating_object=bldg)
    for p in (site, bldg, storey):
        api.geometry.edit_object_placement(f, product=p)
    mat = api.material.add_material(f, name="Concrete C30/37", category="concrete")

    def faceset(pts, tri, closed):
        pl = f.createIfcCartesianPointList3D([tuple(map(float, p)) for p in pts])
        return f.createIfcTriangulatedFaceSet(pl, None, closed, [tuple(int(i) + 1 for i in t) for t in tri], None)

    def product(cls, nm, items, predef, obj_type=None):
        p = api.root.create_entity(f, ifc_class=cls, name=nm, predefined_type=predef)
        if obj_type:
            p.ObjectType = obj_type
        api.geometry.edit_object_placement(f, product=p)
        rep = f.createIfcShapeRepresentation(body, "Body", "Tessellation", items)
        api.geometry.assign_representation(f, product=p, representation=rep)
        api.spatial.assign_container(f, products=[p], relating_structure=storey)
        api.material.assign_material(f, products=[p], material=mat)
        return p

    pts, tri = shell_solid(m, s)
    v_shell = signed_volume(pts, tri)
    slab = product("IfcSlab", f"Shell {name}", [faceset(pts, tri, True)], "ROOF")
    v_rib = 0.0
    for r, (bp, bt) in zip(m.ribs, rib_boxes(m, s)):
        v_rib += signed_volume(bp, bt)
        mem = product("IfcMember", f"Rib {name}-{r.rib_id + 1}", [faceset(bp, bt, None)], "USERDEFINED", "Rib")
        ps = api.pset.add_pset(f, product=mem, name="Pset_MemberCommon")
        api.pset.edit_pset(f, pset=ps, properties={"LoadBearing": True, "IsExternal": True})
        q = api.pset.add_pset(f, product=mem, name="RibShell_Rib")
        api.pset.edit_pset(f, pset=q, properties={"Width_m": r.width, "Depth_m": float(r.depth),
                                                  "Length_m": float(m.rib_length(r.rib_id))})
    vol = m.volume()
    ps = api.pset.add_pset(f, product=slab, name="Pset_SlabCommon")
    api.pset.edit_pset(f, pset=ps, properties={"LoadBearing": True, "IsExternal": True})
    q = api.pset.add_pset(f, product=slab, name="RibShell_Study")
    props = {"Variant": name, "Thickness": meta["thickness"], "ConcreteVolume_m3": float(vol),
             "EmbodiedCarbonA1A3_tCO2e": float(vol * RHO * ECF / 1e3), "CarbonFactor_kgCO2e_per_kg": ECF,
             "CarbonFactorSource": "IStructE 2022 HTCEC 2nd ed. Table 2.3 (ICE 2020), C32/40 25% GGBS",
             "LambdaC_worst_sway": float(meta["lambda_c"]), "WorstSwayDirection_deg": int(meta["worst_angle"]),
             "LateralAcceleration_g": 0.2, "StudyMesh_h_m": h, "Analysis": "CalculiX 2.23 *BUCKLE, S8R + B32R"}
    if meta.get("lambda_g"):
        props["LambdaG_gravity"] = float(meta["lambda_g"])
    if "rib_length_m" in meta:
        props.update({"RibLength_m": float(meta["rib_length_m"]), "RibCrossings": int(meta["crossings"])})
    api.pset.edit_pset(f, pset=q, properties=props)
    EXP.mkdir(exist_ok=True)
    path = EXP / f"P2_{name}.ifc"
    f.write(str(path))
    # check: re-open, schema-validate, compare solid volume with the FE volume
    import ifcopenshell.validate
    g = ifcopenshell.open(str(path))
    log = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(g, log)
    return {"variant": name, "file": path.name, "schema_errors": len(log.statements),
            "V_fe_m3": vol, "V_ifc_m3": v_shell + v_rib, "n_members": len(m.ribs)}


# ---------------------------------------------------------------------------- DXF
def write_dxf(s, lay: dict, meta: dict) -> None:
    import ezdxf
    from matplotlib import pyplot as plt

    doc = ezdxf.new("R2018", setup=True)
    doc.header["$INSUNITS"] = 6                       # metres
    msp = doc.modelspace()
    colours = {"B": 30, "C": 3, "D": 1, "Dd": 6, "Dx": 5}
    doc.layers.add("SHELL_EDGE", color=7)
    doc.layers.add("SHELL_CONTOUR_0.3m", color=8)
    msp.add_lwpolyline([(0, 0), (12, 0), (12, 12), (0, 12)], close=True, dxfattribs={"layer": "SHELL_EDGE"})
    g = np.linspace(0, 12, 121)
    X, Y = np.meshgrid(g, g)
    cs = plt.contour(X, Y, s.z(X.ravel(), Y.ravel()).reshape(X.shape), levels=np.arange(0.3, 2.4, 0.3))
    for lev, segs in zip(cs.levels, cs.allsegs):
        for sg in segs:
            if len(sg) > 1:
                msp.add_polyline3d([(x, y, float(lev)) for x, y in sg], dxfattribs={"layer": "SHELL_CONTOUR_0.3m"})
    plt.close("all")
    for n, c in colours.items():
        if n not in lay:
            continue
        ly = f"RIBS_{n}"
        doc.layers.add(ly, color=c)
        for ln in lay[n]:
            ln = np.asarray(ln)
            msp.add_polyline3d([(float(x), float(y), float(s.z(x, y))) for x, y in ln], dxfattribs={"layer": ly})
        if n in meta:
            msp.add_text(f"{n}: {meta[n]['n_ribs']} ribs, {meta[n]['rib_length_m']:.1f} m, "
                         f"{meta[n]['rib_depth_mm']:.0f} mm deep x 150 mm",
                         dxfattribs={"layer": ly, "height": 0.18}).set_placement((0, -0.6 - 0.35 * list(colours).index(n)))
    if "E" in meta:
        doc.layers.add("E_ZONES", color=4)
        tz = np.asarray(meta["E"]["zone_thickness_m"])
        nz = int(round(np.sqrt(len(tz))))
        for i in range(nz + 1):
            msp.add_line((i * 12 / nz, 0), (i * 12 / nz, 12), dxfattribs={"layer": "E_ZONES"})
            msp.add_line((0, i * 12 / nz), (12, i * 12 / nz), dxfattribs={"layer": "E_ZONES"})
        for ix in range(nz):
            for iy in range(nz):
                msp.add_text(f"{tz[ix * nz + iy] * 1e3:.0f} mm", dxfattribs={"layer": "E_ZONES", "height": 0.25}) \
                    .set_placement(((ix + 0.5) * 12 / nz, (iy + 0.5) * 12 / nz),
                                   align=ezdxf.enums.TextEntityAlignment.MIDDLE_CENTER)
    msp.add_text("Stress-line ribs vs optimised thickness - plan layouts at V = V60 (units: m)",
                 dxfattribs={"layer": "SHELL_EDGE", "height": 0.3}).set_placement((0, 12.6))
    EXP.mkdir(exist_ok=True)
    doc.saveas(EXP / "P2_layouts.dxf")
    ezdxf.readfile(EXP / "P2_layouts.dxf").audit()


# ---------------------------------------------------------------------------- renders
def renders(ms: dict[str, Model], s, modes: dict) -> None:
    import pyvista as pv

    pv.OFF_SCREEN = True
    FIG.mkdir(exist_ok=True)

    def poly(xyz, quads):
        return pv.PolyData(xyz, np.hstack([np.full((len(quads), 1), 4), quads]).ravel())

    def rib_mesh(m):
        blocks = []
        for bp, bt in rib_boxes(m, s):
            blocks.append(pv.PolyData(bp, np.hstack([np.full((len(bt), 1), 3), bt]).ravel()))
        return blocks

    def scene(pl):
        pl.set_background("#fcfcfb")
        pl.enable_anti_aliasing("ssaa")
        pl.camera_position = [(24, -14, 14), (6, 6, 0.8), (0, 0, 1)]
        pl.camera.zoom(1.2)

    for k, md in modes.items():
        pl = pv.Plotter(off_screen=True, window_size=(1600, 1100), lighting="three lights")
        scene(pl)
        sc = 0.6
        mesh = poly(md["xyz"] + sc * md["u"], md["quads"])
        mesh["|u| / max"] = md["amp"]
        pl.add_mesh(mesh, scalars="|u| / max", cmap=MODE_CMAP, smooth_shading=True, show_edges=False,
                    scalar_bar_args={"title": "mode 1 amplitude", "color": "#0b0b0b", "vertical": True,
                                     "position_x": 0.9, "position_y": 0.1, "height": 0.5})
        pl.add_mesh(poly(md["xyz"], md["quads"]), style="wireframe", color="#b9b8b2", opacity=0.25)
        if k in ms:                                 # rib lines drawn on the undeformed surface
            m = ms[k]
            idx = {n: i for i, n in enumerate(sorted({n for c in m.mesh.quads.values() for n in c[:4]}))}
            for r in m.ribs:
                chain = [a for a, _, _ in m.mesh.ribs[r.rib_id]] + [m.mesh.ribs[r.rib_id][-1][2]]
                p = np.array([m.mesh.nodes[n] for n in chain]) + np.array([0, 0, 0.04])
                pl.add_mesh(pv.lines_from_points(p), color="#0b0b0b", line_width=4)
        pl.add_text(f"{k}: mode 1, sway {md['angle']} deg, lambda = {md['factor']:.1f}", color="#0b0b0b",
                    font_size=14)
        pl.screenshot(FIG / f"render_mode_{k}.png")
        pl.close()

    if "D" in ms:                                   # the ribbed shell seen from below
        m = ms["D"]
        pl = pv.Plotter(off_screen=True, window_size=(1600, 1100), lighting="three lights")
        pl.set_background("#fcfcfb")
        pl.enable_anti_aliasing("ssaa")
        pts, tri = shell_solid(m, s)
        pl.add_mesh(pv.PolyData(pts, np.hstack([np.full((len(tri), 1), 3), tri]).ravel()), color="#d9d6cf",
                    smooth_shading=True)
        for b in rib_mesh(m):
            pl.add_mesh(b, color="#eda100")
        pl.camera_position = [(-9, -9, -9), (6, 6, 1.2), (0, 0, 1)]
        pl.camera.zoom(1.25)
        pl.add_text("D: 50 mm shell with sway-envelope ribs, seen from below", color="#0b0b0b", font_size=14)
        pl.screenshot(FIG / "render_D_below.png")
        pl.close()

    if "D" in modes:                                # mode-1 animation
        md = modes["D"]
        pl = pv.Plotter(off_screen=True, window_size=(1000, 720))
        scene(pl)
        mesh = poly(md["xyz"], md["quads"])
        mesh["amp"] = md["amp"]
        pl.add_mesh(mesh, scalars="amp", cmap=MODE_CMAP, clim=(0, 1), smooth_shading=True, show_scalar_bar=False)
        pl.add_text(f"D, buckling mode 1 under sway (lambda = {md['factor']:.1f})", color="#0b0b0b", font_size=12)
        pl.open_gif(str(FIG / "mode_D.gif"), fps=15)
        for ph in np.linspace(0, 2 * np.pi, 36, endpoint=False):
            mesh.points = md["xyz"] + 0.9 * np.sin(ph) * md["u"]
            pl.write_frame()
        pl.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.4)
    ap.add_argument("--skip-modes", action="store_true", help="reuse results/modes_h<h>.npz")
    h = ap.parse_args().h
    a = ap.parse_args()
    s, ms, meta = build_models(h)
    npz = RES / f"modes_h{h:g}.npz"
    if a.skip_modes and npz.exists():
        z = np.load(npz)
        modes = {k: {f: z[f"{k}_{f}"] for f in ("xyz", "quads", "u", "amp", "factor", "angle")}
                 for k in ("A", "C", "D", "E") if f"{k}_xyz" in z}
    else:
        modes = compute_modes(ms, meta, h)
    checks = [write_ifc(k, ms[k], s, meta[k], h) for k in ms]
    lines = ["IFC4 export check: schema errors from ifcopenshell.validate; closed-solid volume vs FE volume",
             "(the FE volume counts rib crossings twice and the rib boxes overlap there in the same way)", ""]
    for c in checks:
        lines.append(f"{c['file']:12s} errors {c['schema_errors']:3d}  members {c['n_members']:3d}  "
                     f"V_fe {c['V_fe_m3']:.4f}  V_ifc {c['V_ifc_m3']:.4f}  diff {100 * (c['V_ifc_m3'] / c['V_fe_m3'] - 1):+.2f} %")
    (EXP / "ifc_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    write_dxf(s, json.load(open(RES / "layouts.json")), meta)
    renders(ms, s, modes)
    print("exports written to", EXP, "and", FIG)


if __name__ == "__main__":
    main()
