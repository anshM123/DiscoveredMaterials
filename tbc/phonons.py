"""Phonon engine: MACE-OMAT-0 force constants -> interface-normal spectral flux.

For each crystal and each interface normal n we compute, over the full
Brillouin zone (no symmetry reduction, so that velocity projections are exact):

    D(w)      = 1/(V N_q) sum_{q,j} delta(w - w_qj)                    [s m^-3]
    M_n(w)    = 1/(2 V N_q) sum_{q,j} |v_qj . n| delta(w - w_qj)       [m^-2]
    Mperp_n(w)= same, weighted by the fraction of the mode displacement
                along n (sum_atoms |e_a . n|^2)                          [m^-2]

M_n(w) dw is the number of phonon channels per unit area and unit time
crossing a plane of normal n from one side (spectral "phonon radiation"
flux), i.e. the quantity that enters the Landauer conductance. For an
isotropic Debye solid it reduces to sum_j D_j v_j / 4.

w is angular frequency (rad/s). Results are cached in data/phonons/<name>.npz.
"""
import os, sys, warnings, json
import numpy as np

# silence third-party library noise only; this package's own warnings
# (missing e-ph data, estimated g, fallbacks) must reach the user
for _cat in (DeprecationWarning, FutureWarning):
    warnings.filterwarnings("ignore", category=_cat)
warnings.filterwarnings("ignore", module=r"(torch|e3nn|mace|ase|phonopy|spglib|numpy|torch_dftd)(\.|$)")
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "data", "phonons")
MODEL = "medium-omat-0"
AMU = 1.66053906660e-27
NW = 600            # frequency bins
W_MAX = 2 * np.pi * 45e12   # rad/s upper edge of grid (45 THz covers diamond)
W_EDGES = np.linspace(0.0, W_MAX, NW + 1)
W = 0.5 * (W_EDGES[1:] + W_EDGES[:-1])
DW = W_EDGES[1] - W_EDGES[0]
T_GRID = np.arange(10.0, 1010.0, 10.0)

_calc = None


def calculator():
    global _calc
    if _calc is None:
        from mace.calculators import mace_mp
        _calc = mace_mp(model=MODEL, default_dtype="float64", device="cpu")
    return _calc


def relax(atoms):
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS
    atoms = atoms.copy()
    atoms.calc = calculator()
    BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=2e-4, steps=500)
    return atoms


def build_phonopy(name):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    from materials import CRYSTALS
    spec = CRYSTALS[name]
    a = relax(spec["build"]())
    u = PhonopyAtoms(symbols=a.get_chemical_symbols(), cell=a.cell[:],
                     scaled_positions=a.get_scaled_positions())
    ph = Phonopy(u, supercell_matrix=np.diag(spec["sc"]), primitive_matrix=None)
    ph.generate_displacements(distance=0.01)
    calc = calculator()
    forces = []
    for sc in ph.supercells_with_displacements:
        s = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        s.calc = calc
        f = s.get_forces()
        forces.append(f - f.mean(axis=0))   # remove residual drift
    ph.forces = forces
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    if "nac" in spec:
        # isotropic Born charges: cation +Z, anions balance the charge
        Z, eps = spec["nac"]["Z"], spec["nac"]["eps"]
        syms = list(ph.primitive.symbols)
        cation = spec["nac"].get("cation", syms[0])
        ncat = sum(s == cation for s in syms); nan = len(syms) - ncat
        born = np.array([np.eye(3) * (Z if s == cation else -Z * ncat / nan) for s in syms])
        # factor = e^2/(4 pi eps0) in eV*Angstrom (Hartree*Bohr), for eV/Angstrom force units
        ph.nac_params = {"born": born, "dielectric": np.eye(3) * eps, "factor": 14.399645}
    return ph, a


def sound_velocities(ph, ndir=200, dq=0.01, seed=0):
    """Debye-averaged LA and TA velocities from acoustic group velocities at small q."""
    rng = np.random.default_rng(seed)
    rec = np.linalg.inv(ph.primitive.cell).T          # rows = b_i / (2 pi), 1/Angstrom
    dirs = rng.normal(size=(ndir, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    # q (reduced) such that cartesian q = dq * dir (1/Angstrom, without 2 pi)
    qred = np.array([np.linalg.solve(rec.T, dq * d) for d in dirs])
    ph.run_qpoints(qred, with_group_velocities=True)
    gv = ph.get_qpoints_dict()["group_velocities"]    # (nq, nband, 3) THz*Angstrom
    speeds = np.linalg.norm(gv[:, :3, :], axis=2) * 100.0   # m/s
    speeds.sort(axis=1)
    vT = np.mean(speeds[:, :2] ** -3) ** (-1 / 3)
    vL = np.mean(speeds[:, 2] ** -3) ** (-1 / 3)
    return float(vL), float(vT)


def compute(name, force=False):
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, f"{name}.npz")
    if os.path.exists(out) and not force:
        return load(name)
    from materials import CRYSTALS
    spec = CRYSTALS[name]
    ph, a = build_phonopy(name)
    V = abs(np.linalg.det(ph.primitive.cell)) * 1e-30      # m^3 per primitive cell
    nat = len(ph.primitive.masses)
    mass = float(np.sum(ph.primitive.masses)) * AMU
    ph.run_mesh(spec["mesh"], with_group_velocities=True, with_eigenvectors=True,
                is_mesh_symmetry=False, is_gamma_center=True)
    md = ph.get_mesh_dict()
    f = md["frequencies"]                                   # THz (nq, nb)
    gv = md["group_velocities"] * 100.0                     # m/s (nq, nb, 3)
    ev = md["eigenvectors"]                                 # (nq, 3*nat, nb)
    nq = f.shape[0]
    w = 2 * np.pi * np.clip(f, 0, None) * 1e12
    idx = np.clip((w / DW).astype(int), 0, NW - 1)
    D = np.bincount(idx.ravel(), minlength=NW) / (V * nq * DW)
    res = dict(W=W, D=D, V=V, nat=nat, rho=mass / V, n_atoms=nat / V,
               m_avg=mass / nat / AMU, a_cell=np.array(a.cell[:]),
               f_max=float(f.max()))
    ev_r = ev.reshape(nq, nat, 3, -1)                       # (nq, atom, xyz, band)
    for lab, nvec in spec["normals"].items():
        n = np.asarray(nvec, float)
        n /= np.linalg.norm(n)
        vn = np.abs(gv @ n)                                 # (nq, nb)
        M = np.bincount(idx.ravel(), weights=vn.ravel(), minlength=NW) / (2 * V * nq * DW)
        proj = np.abs(np.einsum("qaxb,x->qab", ev_r, n)) ** 2
        wperp = proj.sum(axis=1)                            # (nq, nb)
        Mp = np.bincount(idx.ravel(), weights=(vn * wperp).ravel(), minlength=NW) / (2 * V * nq * DW)
        res[f"M_{lab}"] = M
        res[f"Mperp_{lab}"] = Mp
    del md, ev, ev_r
    # Mean-square displacement vs T (Angstrom^2, averaged over atoms, full 3D)
    ph.run_mesh([max(8, m // 2) for m in spec["mesh"]], with_eigenvectors=True,
                is_mesh_symmetry=False)
    ph.run_thermal_displacements(temperatures=T_GRID)
    td = ph.get_thermal_displacements_dict()["thermal_displacements"]  # (nT, nat, 3)
    res["T"] = T_GRID
    res["msd"] = np.asarray(td).reshape(len(T_GRID), -1).sum(axis=1) / nat
    res["vL"], res["vT"] = sound_velocities(ph)
    np.savez(out, **res)
    ph.save(os.path.join(CACHE, f"{name}_phonopy.yaml"), settings={"force_constants": True})
    return load(name)


def load(name):
    d = np.load(os.path.join(CACHE, f"{name}.npz"))
    return {k: d[k] for k in d.files}


if __name__ == "__main__":
    from materials import CRYSTALS
    names = sys.argv[1:] or list(CRYSTALS)
    for nm in names:
        r = compute(nm, force=True)
        Ms = {k: float(np.sum(r[k]) * DW) for k in r if k.startswith("M_")}
        print(json.dumps(dict(name=nm, f_max_THz=round(float(r["f_max"]), 2),
                              vL=round(float(r["vL"])), vT=round(float(r["vT"])),
                              rho=round(float(r["rho"])),
                              msd300=round(float(np.interp(300, r["T"], r["msd"])), 5),
                              dos_norm=round(float(np.sum(r["D"]) * DW * r["V"] / (3 * r["nat"])), 4),
                              **{k: f"{v:.3e}" for k, v in Ms.items()})), flush=True)


# --------------------------------------------------------------------------
# k_par-resolved channel flux for specular (momentum-conserving) transmission
CH_CACHE = os.path.join(CACHE, "channels")


def channel_spectra(name, normal, nk=24, nz=24, force=False, with_velocity=False):
    """M_k(w): the spectral flux of phonons crossing a plane of normal n,
    resolved by in-plane wavevector k_par on an nk x nk grid of the reduced
    2D reciprocal cell of the interface-plane Bravais lattice. Evaluated on
    columns q = k_par + q_z n, q_z in [0, 2 pi/d), d = V_prim / A_plane, with
    the cached MACE force constants (any orientation). sum_k M_k(w) / nk^2 =
    M_n(w) (checked). Returns array (nk*nk, NW) in m^-2 per (rad/s)."""
    from phonopy import load as phload
    from materials import CRYSTALS
    from model import plane_lattice
    os.makedirs(CH_CACHE, exist_ok=True)
    out = os.path.join(CH_CACHE, f"{name}_{normal}_{nk}_{nz}.npz")
    if os.path.exists(out) and not force:
        z = np.load(out)
        return (z["Mk"], z["Vk"]) if with_velocity else z["Mk"]
    ph = phload(os.path.join(CACHE, f"{name}_phonopy.yaml"), produce_fc=False, log_level=0)
    spec = CRYSTALS[name]
    if "nac" in spec:
        pass   # NAC already stored in the yaml if it was set
    d = load(name)
    n = np.asarray(spec["normals"][normal], float); n /= np.linalg.norm(n)
    b2 = plane_lattice(d["a_cell"], n)                       # 2x2 in plane frame
    cellp = np.asarray(ph.primitive.cell)
    V = abs(np.linalg.det(cellp))
    A = abs(np.linalg.det(b2))
    dperp = V / A
    # plane frame -> Cartesian
    # recover e1,e2 used by plane_lattice
    cell = np.asarray(d["a_cell"])
    r = np.arange(-4, 5)
    I3 = np.array(np.meshgrid(r, r, r)).reshape(3, -1).T
    Vv = I3 @ cell; L = np.linalg.norm(Vv, axis=1)
    ok = (L > 1e-6) & (np.abs(Vv @ n) < 1e-4 * np.maximum(L, 1))
    v1 = Vv[ok][np.argmin(L[ok])]
    e1 = v1 / np.linalg.norm(v1); e2 = np.cross(n, e1)
    rec2 = 2 * np.pi * np.linalg.inv(b2).T                   # rows: 2D reciprocal vectors (plane frame)
    f = (np.arange(nk) + 0.5) / nk
    K2 = np.array([[a, b] for a in f for b in f]) @ rec2     # plane-frame k_par
    Kc = K2[:, :1] * e1 + K2[:, 1:2] * e2                    # Cartesian, 1/Angstrom (with 2 pi)
    qz = (np.arange(nz) + 0.5) / nz * 2 * np.pi / dperp
    Q = (Kc[:, None, :] + qz[None, :, None] * n).reshape(-1, 3)
    qred = Q @ cellp.T / (2 * np.pi)                          # reduced coords of primitive reciprocal
    ph.run_qpoints(qred)
    fr = ph.get_qpoints_dict()["frequencies"].reshape(nk * nk, nz, -1)
    w = 2 * np.pi * np.clip(fr, 0, None) * 1e12                    # (nk2, nz, nb)
    # exact channel counting by level crossings: along each periodic q_z
    # column, the number of right-moving channels at frequency w equals half
    # the number of times the (piecewise-linear) band crosses w.
    w2 = np.concatenate([w, w[:, :1]], axis=1)
    lo = np.minimum(w2[:, :-1], w2[:, 1:]); hi = np.maximum(w2[:, :-1], w2[:, 1:])
    ilo = np.clip(np.round(lo / DW).astype(int), 0, NW)
    ihi = np.clip(np.round(hi / DW).astype(int), 0, NW)
    diff = np.zeros((nk * nk, NW + 1))
    kk = np.broadcast_to(np.arange(nk * nk)[:, None, None], ilo.shape)
    np.add.at(diff, (kk.ravel(), ilo.ravel()), 1.0)
    np.add.at(diff, (kk.ravel(), ihi.ravel()), -1.0)
    Nk = np.cumsum(diff, axis=1)[:, :NW] / 2.0                     # right-moving channels
    # mean |v_z| of the crossing modes at (k, w): segment slope dw/dq_z
    dq = 2 * np.pi / dperp / nz * 1e10                             # 1/m
    vseg = (hi - lo) / dq                                          # m/s
    dv = np.zeros((nk * nk, NW + 1))
    np.add.at(dv, (kk.ravel(), ilo.ravel()), vseg.ravel())
    np.add.at(dv, (kk.ravel(), ihi.ravel()), -vseg.ravel())
    Vk = np.divide(np.cumsum(dv, axis=1)[:, :NW] / 2.0, Nk, out=np.zeros_like(Nk), where=Nk > 0)
    # Landauer: M(w) = (1/(2 pi A)) <N_k(w)>_k  ->  per-k contribution
    Mk = Nk / (2 * np.pi * A * 1e-20)
    np.savez_compressed(out, Mk=Mk.astype(np.float32), Vk=Vk.astype(np.float32))
    return (Mk, Vk) if with_velocity else Mk

