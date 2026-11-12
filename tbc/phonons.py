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

warnings.filterwarnings("ignore")
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
        Z, eps = spec["nac"]["Z"], spec["nac"]["eps"]
        syms = a.get_chemical_symbols()
        cation = syms[0]
        born = np.array([np.eye(3) * (Z if s == cation else -Z) for s in syms])
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
