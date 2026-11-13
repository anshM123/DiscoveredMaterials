"""Oxide glasses from a MACE-OMAT-0 melt-quench, and their propagating
(propagon) phonon flux for the fast tier (model.Side, GLASS_FLUX).

Pipeline (python glass.py <name>), cached in data/glass/:
  1. make_glass: random atoms at the experimental density in a hexagonal cell
     (a = A_LAT, laterally matched to Al(111) 4x4 and Si(111) 3x3 for the AGF
     references of agf_glass.py) -> Langevin MD at T_melt (t_melt) -> linear
     quench to 300 K (t_quench, ~1e15 K/s, typical of MD glasses) -> 300 K
     (t_anneal) -> FIRE relaxation (float32, then float64).
     MACE-OMAT-0 without D3 (covalent/ionic oxide).
  2. force_constants: phonopy finite displacements (0.01 A, float64 MACE)
     of the periodic glass -> <prefix>_ph.yaml.
  3. flux: propagating spectral flux M(w) = 1/(2V) sum |v_g| delta(w - w_q)
     on a 6x6x6 q-mesh (average of the normal and in-plane directions; the
     glass is isotropic) -> <prefix>_flux.npz. Diffusons (flat bands of the
     periodic glass) carry ~no flux here; model.Side adds them back with the
     global weight W_DIFFUSON.
  Also reported: coordination and relaxed-ion sound velocities.

The cached a-SiO2 (aSiO2_32) used T_melt = 4000 K, t_melt = 2 ps,
t_quench = 4 ps, t_anneal = 0.5 ps, seed 0 (the defaults below).
"""
import os, sys, json, warnings
import numpy as np

# silence third-party library noise only; this package's own warnings
# (missing e-ph data, estimated g, fallbacks) must reach the user
for _cat in (DeprecationWarning, FutureWarning):
    warnings.filterwarnings("ignore", category=_cat)
warnings.filterwarnings("ignore", module=r"(torch|e3nn|mace|ase|phonopy|spglib|numpy|torch_dftd)(\.|$)")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "glass")
A_LAT = 11.45            # Angstrom: Al(111) 4x4 = 11.41, Si(111) 3x3 = 11.51

# name (as in materials.AMORPHOUS) -> composition, molar mass (g/mol per
# formula unit), density of the amorphous solid (g/cm^3), number of formula units
GLASSES = {
    "SiO2": dict(prefix="aSiO2", formula=(("Si", 1), ("O", 2)), M=60.084, rho=2.20, units=32,
                 dmin={"Si": 1.6, "O": 1.6}),
    # amorphous (ALD / anodic) alumina: density 3.0-3.1 g/cm^3 (ALD films)
    "aAl2O3": dict(prefix="aAl2O3", formula=(("Al", 2), ("O", 3)), M=101.961, rho=3.05, units=20,
                   dmin={"Al": 1.6, "O": 1.6}),
}


def _calcs():
    from mace.calculators import mace_mp
    return (lambda: mace_mp(model="medium-omat-0", default_dtype="float32", device="cpu"),
            lambda: mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu"))


def make_glass(name="SiO2", n_units=None, seed=0, T_melt=4000.0, t_melt=2.0, t_quench=4.0, t_anneal=0.5, dt=1.0):
    from ase import Atoms, units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.optimize import FIRE
    g = GLASSES[name]
    n_units = n_units or g["units"]
    calc32, calc = _calcs()
    rng = np.random.default_rng(seed)
    mass = n_units * g["M"] / 6.02214076e23            # g
    vol = mass / g["rho"] * 1e24                       # A^3
    area = A_LAT ** 2 * np.sqrt(3) / 2
    c = vol / area
    cell = [[A_LAT, 0, 0], [-A_LAT / 2, A_LAT * np.sqrt(3) / 2, 0], [0, 0, c]]
    # random placement with a minimum distance
    pos, sym = [], []
    for s_, n in g["formula"]:
        for _ in range(n * n_units):
            for _try in range(2000):
                p = rng.random(3) @ np.array(cell)
                ok = True
                if pos:
                    a = Atoms(sym + [s_], positions=pos + [p], cell=cell, pbc=True)
                    d = a.get_distances(len(a) - 1, list(range(len(a) - 1)), mic=True)
                    ok = d.min() > g["dmin"][s_]
                if ok:
                    break
            pos.append(p); sym.append(s_)
    at = Atoms(sym, positions=pos, cell=cell, pbc=True)
    at.calc = calc32()
    MaxwellBoltzmannDistribution(at, temperature_K=T_melt, rng=np.random.default_rng(seed + 1))
    dyn = Langevin(at, dt * units.fs, temperature_K=T_melt, friction=0.02, rng=np.random.default_rng(seed + 2))
    log = lambda msg: print(msg, flush=True)
    n_melt = int(t_melt * 1000 / dt); n_q = int(t_quench * 1000 / dt); n_a = int(t_anneal * 1000 / dt)
    log(f"glass: {len(at)} atoms, c = {c:.2f} A; melt {n_melt} steps")
    dyn.run(n_melt)
    nseg = 24
    for k in range(nseg):
        T = T_melt + (300.0 - T_melt) * (k + 1) / nseg
        dyn.set_temperature(temperature_K=T)
        dyn.run(n_q // nseg)
        if k % 6 == 5:
            log(f"  quench T = {T:.0f} K")
    dyn.set_temperature(temperature_K=300.0)
    dyn.run(n_a)
    at.set_velocities(np.zeros((len(at), 3)))
    FIRE(at, logfile=None).run(fmax=0.05, steps=1000)
    at.calc = calc()
    FIRE(at, logfile=None).run(fmax=0.01, steps=1000)
    return at


def analyse(at, cutoff=None):
    """Mean coordination of each cation and of O (bond cutoff 2.0 A for Si-O,
    2.3 A otherwise), fraction of the dominant coordination, and density."""
    from ase.neighborlist import neighbor_list
    sy = np.array(at.get_chemical_symbols())
    cutoff = cutoff or (2.0 if "Si" in sy else 2.3)
    i, j = neighbor_list("ij", at, cutoff)
    cn = np.bincount(i, minlength=len(at))
    out = dict(density=float(np.sum(at.get_masses()) / 6.02214076e23 / at.get_volume() * 1e24))
    for el in sorted(set(sy)):
        c = cn[sy == el]
        mode = int(np.bincount(c).argmax())
        out[f"{el}_cn_mean"] = float(c.mean())
        out[f"{el}{mode}"] = float(np.mean(c == mode))
    return out


def force_constants(at):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    calc = _calcs()[1]()
    u = PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.cell[:], scaled_positions=at.get_scaled_positions())
    ph = Phonopy(u, supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")
    ph.generate_displacements(distance=0.01)
    F = []
    for sc in ph.supercells_with_displacements:
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True); a.calc = calc
        f = a.get_forces(); F.append(f - f.mean(0))
    ph.forces = F; ph.produce_force_constants(); ph.symmetrize_force_constants()
    return ph


def sound_velocities(at, eps=0.005):
    """Isotropic sound velocities (m/s) from the relaxed-ion elastic constants
    of the glass: +-eps normal and shear strains, ions relaxed (float64 MACE),
    C11 and C44 averaged over the three axes; vL = sqrt(C11/rho), vT = sqrt(C44/rho)."""
    from ase.optimize import FIRE
    calc = _calcs()[1]()
    rho = np.sum(at.get_masses()) * 1.66053906660e-27 / (at.get_volume() * 1e-30)
    cell0 = at.cell.array.copy()

    def stress(e):
        a = at.copy(); a.calc = calc
        a.set_cell(cell0 @ (np.eye(3) + e), scale_atoms=True)
        FIRE(a, logfile=None).run(fmax=0.005, steps=500)
        return a.get_stress(voigt=False) * 160.21766208e9          # eV/A^3 -> Pa

    C11, C44 = [], []
    for i in range(3):
        e = np.zeros((3, 3)); e[i, i] = eps
        C11.append((stress(e)[i, i] - stress(-e)[i, i]) / (2 * eps))
        j = (i + 1) % 3
        e = np.zeros((3, 3)); e[i, j] = e[j, i] = eps / 2               # engineering shear strain eps
        C44.append((stress(e)[i, j] - stress(-e)[i, j]) / (2 * eps))
    return float(np.sqrt(np.mean(C11) / rho)), float(np.sqrt(np.mean(C44) / rho))


def flux(ph, mesh=(6, 6, 6)):
    """Propagating spectral flux of the periodic glass on the phonons.W grid."""
    from phonons import W, DW, NW
    V = abs(np.linalg.det(ph.primitive.cell)) * 1e-30
    ph.run_mesh(list(mesh), with_group_velocities=True, is_mesh_symmetry=False, is_gamma_center=False)
    md = ph.get_mesh_dict(); f = md["frequencies"]; gv = md["group_velocities"] * 100
    w = 2 * np.pi * np.clip(f, 0, None) * 1e12
    idx = np.clip((w / DW).astype(int), 0, NW - 1); nq = f.shape[0]
    D = np.bincount(idx.ravel(), minlength=NW) / (V * nq * DW)
    Mz = np.bincount(idx.ravel(), weights=np.abs(gv[:, :, 2]).ravel(), minlength=NW) / (2 * V * nq * DW)
    Mx = np.bincount(idx.ravel(), weights=np.abs(gv[:, :, 0]).ravel(), minlength=NW) / (2 * V * nq * DW)
    return dict(W=W, M=0.5 * (Mz + Mx), D=D, Mz=Mz, Mx=Mx)


def build(name, n_units=None, seed=0):
    """Run (or resume) the full pipeline for one glass; returns the summary."""
    from ase.io import read, write
    from phonopy import load as phload
    os.makedirs(OUT, exist_ok=True)
    g = GLASSES[name]
    tag = f"{g['prefix']}_{n_units or g['units']}"
    xyz, phf, fx = (os.path.join(OUT, tag + s) for s in (".extxyz", "_ph.yaml", "_flux.npz"))
    if os.path.exists(xyz):
        at = read(xyz)
    else:
        at = make_glass(name, n_units, seed)
        write(xyz, at)
    st = analyse(at)
    print(json.dumps(st), flush=True)
    if os.path.exists(phf):
        ph = phload(phf, log_level=0)
    else:
        ph = force_constants(at)
        ph.save(phf, settings={"force_constants": True})
    if not os.path.exists(fx):
        np.savez(fx, **flux(ph))
    vL, vT = sound_velocities(at)
    st.update(vL_MACE=vL, vT_MACE=vT)
    json.dump(st, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    print(json.dumps(st), flush=True)
    return st


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "SiO2")
