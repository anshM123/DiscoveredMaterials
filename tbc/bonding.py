"""Interfacial bond stiffness from MACE-OMAT-0 + D3(BJ) dispersion.

For an interface we compute the energy of a rigid displacement delta of the
upper part against the lower part, E(delta), around the relaxed contact, and
take the stiffness per unit area

    K_n = E''(delta_z) / A,   K_t = E''(delta_x) / A      [N m^-3]

The reference ("bulk-like bond") stiffness of each material across a lattice
plane is K_bulk = rho v^2 / d with d = n^(-1/3) (v = v_L for K_n, v_T for
K_t) and the interface reference is the harmonic mean of both sides. The
bond factor used in model.Interface is spring(K)/spring(K_ref) <= 1, so a
chemically bonded interface (K >= K_ref) is not penalised.

Results are cached in data/bonding.json.
"""
import json, os, sys, warnings
import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "data", "bonding.json")
EV = 1.602176634e-19
PMMA_UNITS_PER_A2 = (1180 / (100.12e-3 / 6.02214076e23)) ** (2 / 3) * 1e-20

_calc = None
_calc32 = None


def calc():
    global _calc
    if _calc is None:
        from mace.calculators import mace_mp
        _calc = mace_mp(model="medium-omat-0", dispersion=True, default_dtype="float64", device="cpu")
    return _calc


def calc32():
    """float32 MACE + D3: 2.4x faster, used for MD and pre-relaxation only;
    energies for stiffness scans always use float64 (calc())."""
    global _calc32
    if _calc32 is None:
        from mace.calculators import mace_mp
        _calc32 = mace_mp(model="medium-omat-0", dispersion=True, default_dtype="float32", device="cpu")
    return _calc32


REACTIVE = dict(T=700.0, ps=0.5, dt_fs=1.0)


def _reactive_anneal(sys_, nsub, log=None):
    """Reactive protocol for molecule/surface contacts: Langevin MD at
    REACTIVE['T'] for REACTIVE['ps'] (H given deuterium mass so that a 1 fs
    step is stable; masses do not affect static stiffness), then quench:
    float32 relaxation and float64 polish. Lets chemisorption, bond
    breaking and oxidation (e.g. Ni-O, Ti-O/Ti-C formation) happen if they are
    activated at this temperature on a ps time scale. The bottom 60% of the
    substrate is held fixed throughout."""
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS, FIRE
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase import units
    z = sys_.positions[:, 2]
    zs = z[:nsub]
    fixed = [i for i in range(nsub) if zs[i] < zs.min() + 0.6 * (zs.max() - zs.min())]
    # The PMMA repeat units stand for a polymer film: the rest of the film keeps
    # them from desorbing. A one-sided harmonic wall 1.5 A above the initial top
    # of the molecular layer (Hookean plane restraint, k = 5 eV/A^2) represents it.
    from ase.constraints import Hookean
    zwall = sys_.positions[nsub:, 2].max() + 1.5
    walls = [Hookean(a1=i, a2=(0.0, 0.0, 1.0, -zwall), k=5.0) for i in range(nsub, len(sys_))]   # keeps z <= zwall
    sys_.set_constraint([FixAtoms(indices=fixed)] + walls)
    m = sys_.get_masses()
    m[np.array(sys_.get_chemical_symbols()) == "H"] = 2.014
    sys_.set_masses(m)
    sys_.calc = calc32()
    MaxwellBoltzmannDistribution(sys_, temperature_K=REACTIVE["T"], rng=np.random.default_rng(0))
    nsteps = int(REACTIVE["ps"] * 1000 / REACTIVE["dt_fs"])
    dyn = Langevin(sys_, REACTIVE["dt_fs"] * units.fs, temperature_K=REACTIVE["T"], friction=0.02,
                   rng=np.random.default_rng(1))
    dyn.run(nsteps)
    sys_.set_velocities(np.zeros((len(sys_), 3)))
    FIRE(sys_, logfile=None).run(fmax=0.05, steps=600)
    sys_.calc = calc()
    BFGS(sys_, logfile=None).run(fmax=0.05, steps=150)
    sys_.set_constraint()
    return sys_


def relaxed_a(name):
    from phonons import load
    cell = load(name)["a_cell"]
    return cell


def fcc_a(name):
    cell = relaxed_a(name)
    return np.linalg.norm(cell[0]) * np.sqrt(2)


def metal_slab(name, size=(4, 4, 4)):
    from ase.build import fcc111, hcp0001
    if name == "Ti":
        cell = relaxed_a("Ti")
        a, c = np.linalg.norm(cell[0]), np.linalg.norm(cell[2])
        s = hcp0001("Ti", size=size, a=a, c=c, vacuum=0.0, orthogonal=False)
    else:
        s = fcc111(name, size=size, a=fcc_a(name), vacuum=0.0, orthogonal=False)
    return s


def pmma_unit(seed):
    """Methyl isobutyrate = PMMA repeat unit with H caps, random orientation."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    m = Chem.AddHs(Chem.MolFromSmiles("COC(=O)C(C)C"))
    AllChem.EmbedMolecule(m, randomSeed=int(seed) + 1)
    AllChem.MMFFOptimizeMolecule(m)
    a = Atoms([x.GetSymbol() for x in m.GetAtoms()], positions=m.GetConformer().GetPositions())
    a.positions -= a.positions.mean(axis=0)
    rng = np.random.default_rng(seed)
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, z = q
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    a.positions = a.positions @ R.T
    return a


def hydroxylated_silica(size=(2, 2), nlayers=3):
    """Hydroxylated beta-cristobalite(111) slab (~4.5 OH/nm^2, the silanol
    density of amorphous silica surfaces). Bottom also hydroxylated and fixed."""
    from ase.build import diamond111
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    si = diamond111("Si", size=(size[0], size[1], 2 * nlayers + 2), a=7.16, vacuum=0.0, orthogonal=False)
    si.pbc = (True, True, False)
    si.cell[2] = [0, 0, 60.0]
    # drop the outermost atoms with three dangling bonds -> shuffle-plane
    # termination, one dangling bond per surface Si (capped by OH)
    i, j = neighbor_list("ij", si, 3.2)
    si = si[np.bincount(i, minlength=len(si)) > 1]
    i, j, D = neighbor_list("ijD", si, 3.2)
    pos = list(si.positions); sym = list(si.get_chemical_symbols())
    for a_, b_, d in zip(i, j, D):
        if a_ < b_:
            pos.append(si.positions[a_] + d / 2); sym.append("O")
    nn = np.bincount(i, minlength=len(si))
    zmid = si.positions[:, 2].mean()
    for k in range(len(si)):
        if nn[k] < 4:
            sgn = 1.0 if si.positions[k, 2] > zmid else -1.0
            o = si.positions[k] + np.array([0, 0, sgn * 1.62])
            pos.append(o); sym.append("O")
            pos.append(o + np.array([0.9, 0, sgn * 0.35])); sym.append("H")
    s = Atoms(sym, positions=pos, cell=si.cell, pbc=(True, True, False))
    return s


def _scan(sys_atoms, upper_mask, A, h=0.04):
    """Second derivatives of E w.r.t. rigid shift of upper part (z and x)."""
    c = calc()
    base = sys_atoms.positions.copy()

    def E(dv):
        s = sys_atoms.copy(); s.calc = c
        p = base.copy(); p[upper_mask] += dv
        s.positions = p
        return s.get_potential_energy()

    e0 = E(np.zeros(3))
    out = {}
    for lab, ax in [("n", np.array([0, 0, 1.0])), ("tx", np.array([1.0, 0, 0])), ("ty", np.array([0, 1.0, 0]))]:
        ep, em = E(h * ax), E(-h * ax)
        out[lab] = (ep + em - 2 * e0) / h ** 2 * EV / 1e-20 / (A * 1e-20)   # J/m^2 / m^2 = N/m^3
    return dict(K_n=out["n"], K_t=0.5 * (out["tx"] + out["ty"]))


def _place_and_relax(sub, mols, z_gap=2.4, fmax=0.05, steps=400):
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    sys_ = sub.copy()
    top = sub.positions[:, 2].max()
    nsub = len(sub)
    cell = sub.cell
    n = len(mols)
    nx = int(np.ceil(np.sqrt(n)))
    for k, m in enumerate(mols):
        fx, fy = ((k % nx) + 0.5) / nx, ((k // nx) + 0.5) / nx
        xy = fx * cell[0] + fy * cell[1]
        m = m.copy()
        m.positions += np.array([xy[0], xy[1], top + z_gap - m.positions[:, 2].min()])
        sys_ += m
    sys_.cell[2] = [0, 0, top + 30.0]
    sys_.pbc = (True, True, False)
    zsub = sub.positions[:, 2]
    fixed = [i for i in range(nsub) if zsub[i] < zsub.min() + 0.6 * (zsub.max() - zsub.min())]
    sys_.set_constraint(FixAtoms(indices=fixed))
    sys_.calc = calc()
    BFGS(sys_, logfile=None).run(fmax=fmax, steps=steps)
    sys_.set_constraint()
    upper = np.zeros(len(sys_), bool); upper[nsub:] = True
    return sys_, upper


def _min_dist(s, nsub, sym):
    from ase.geometry import get_distances
    sy = np.array(s.get_chemical_symbols())
    idx = [i for i in range(nsub, len(s)) if sy[i] == sym]
    sub_idx = [i for i in range(nsub) if sy[i] not in ("O", "H")]
    if not idx or not sub_idx:
        return float("nan")
    D = get_distances(s.positions[sub_idx], s.positions[idx], cell=s.cell, pbc=s.pbc)[1]
    return float(D.min())


def _place_and_relax_fast(sub, mols):
    """Same placement as _place_and_relax, pre-relaxed with float32."""
    global _calc
    keep = _calc
    _calc = calc32()
    try:
        s, up = _place_and_relax(sub, mols, steps=200)
    finally:
        _calc = keep
    return s, up


def K_molecular_layer(substrate, seeds=(0, 1, 2), reactive=False):
    """PMMA repeat units at PMMA areal density on a substrate slab.
    reactive=True adds the anneal of _reactive_anneal before the scan."""
    sub = metal_slab(substrate) if substrate != "SiO2" else hydroxylated_silica()
    sub.pbc = (True, True, False)
    A = np.linalg.norm(np.cross(sub.cell[0], sub.cell[1]))
    nmol = max(1, int(round(PMMA_UNITS_PER_A2 * A)))
    res = []
    for sd in seeds:
        mols = [pmma_unit(1000 * sd + k) for k in range(nmol)]
        if reactive:
            s, up = _place_and_relax_fast(sub, mols)
            s = _reactive_anneal(s, len(sub))
            from collections import Counter
            print(f"    {substrate}/PMMA seed {sd} after anneal: min metal/Si-O dist "
                  f"{_min_dist(s, len(sub), 'O'):.2f} A, min substrate-C dist {_min_dist(s, len(sub), 'C'):.2f} A", flush=True)
        else:
            s, up = _place_and_relax(sub, mols)
        r = _scan(s, up, A)
        res.append(r)
        print(f"  {substrate}/PMMA seed {sd}: n_mol={nmol} K_n={r['K_n']:.3e} K_t={r['K_t']:.3e}", flush=True)
    return dict(K_n=float(np.mean([r["K_n"] for r in res])), K_t=float(np.mean([r["K_t"] for r in res])),
                K_n_all=[r["K_n"] for r in res], K_t_all=[r["K_t"] for r in res], A=A, n_mol=nmol)


def _relax_gap(lower, upper_atoms, A, d0):
    """Find equilibrium separation of two rigid slabs by scanning the gap."""
    c = calc()
    best = None
    for d in np.arange(d0 - 0.6, d0 + 0.9, 0.1):
        s = lower.copy()
        u = upper_atoms.copy()
        u.positions[:, 2] += lower.positions[:, 2].max() + d - u.positions[:, 2].min()
        s += u
        s.calc = c
        e = s.get_potential_energy()
        if best is None or e < best[0]:
            best = (e, d, s)
    return best


def K_coherent(lower, upper, d0, relax_fmax=0.05):
    """Two crystalline slabs with the same lateral cell (upper strained to
    lower beforehand). Gap optimised, then interface atoms relaxed with the
    outer halves fixed, then rigid-shift second derivatives."""
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    lower = lower.copy(); upper = upper.copy()
    for s in (lower, upper):
        s.pbc = (True, True, False)
    A = np.linalg.norm(np.cross(lower.cell[0], lower.cell[1]))
    e, d, s = _relax_gap(lower, upper, A, d0)
    s.cell[2] = [0, 0, s.positions[:, 2].max() + 20]
    nl = len(lower)
    z = s.positions[:, 2]
    zl, zu = lower.positions[:, 2], upper.positions[:, 2]
    fix = [i for i in range(nl) if zl[i] < zl.min() + 0.5 * (zl.max() - zl.min())]
    zu_now = z[nl:]
    fix += [nl + i for i in range(len(upper)) if zu_now[i] > zu_now.min() + 0.5 * (zu_now.max() - zu_now.min())]
    s.set_constraint(FixAtoms(indices=fix))
    s.calc = calc()
    BFGS(s, logfile=None).run(fmax=relax_fmax, steps=300)
    s.set_constraint()
    up = np.zeros(len(s), bool); up[nl:] = True
    r = _scan(s, up, A)
    r.update(gap=float(d), A=float(A))
    return r


def _angle(c):
    return np.degrees(np.arccos(c[0] @ c[1] / np.linalg.norm(c[0]) / np.linalg.norm(c[1])))


def match(slab, target_cell):
    """Strain slab laterally onto target lateral cell (same supercell shape).
    A 60-degree hexagonal basis is first re-expressed as 120 degrees (or vice
    versa) so that the lattice, not just the fractional coordinates, maps."""
    from ase.build import make_supercell
    s = slab.copy()
    a_s, a_t = _angle(s.cell), _angle(np.asarray(target_cell))
    if abs(a_s - 60) < 1 and abs(a_t - 120) < 1:
        s = make_supercell(s, np.array([[1, 0, 0], [-1, 1, 0], [0, 0, 1]]), wrap=True)
    elif abs(a_s - 120) < 1 and abs(a_t - 60) < 1:
        s = make_supercell(s, np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]]), wrap=True)
    slab = s; s = slab.copy()
    cell = s.cell.array.copy()
    cell[0], cell[1] = target_cell[0], target_cell[1]
    s.set_cell(cell, scale_atoms=True)
    return s, float(np.linalg.norm(target_cell[0]) / np.linalg.norm(slab.cell[0]) - 1)


def crystal_slab(name, size, nlayers):
    from ase.build import fcc111, diamond111, surface, make_supercell
    from materials import CRYSTALS
    if name in ("Al", "Cu", "Au", "Ni"):
        return fcc111(name, size=(size[0], size[1], nlayers), a=fcc_a(name), vacuum=0.0, orthogonal=False)
    if name in ("Si", "C"):
        # bulk-terminated (111): one dangling bond per surface atom
        from ase.neighborlist import neighbor_list
        s = diamond111(name, size=(size[0], size[1], nlayers + 2), a=fcc_a(name), vacuum=0.0, orthogonal=False)
        s.pbc = (True, True, False)
        i = neighbor_list("i", s, 1.2 * fcc_a(name) * np.sqrt(3) / 4)
        s = s[np.bincount(i, minlength=len(s)) > 1]
        s.positions[:, 2] -= s.positions[:, 2].min()
        return s
    if name == "TiC":
        from ase.build import bulk
        b = bulk("TiC", "rocksalt", a=fcc_a("TiC"))
        s = surface(b, (1, 1, 1), nlayers, vacuum=0.0)
        s = make_supercell(s, np.diag([size[0], size[1], 1]))
        # Ti-terminated on both faces (the stable TiC(111) termination)
        z = s.positions[:, 2]
        top = z > z.max() - 0.3
        if s.get_chemical_symbols()[int(np.argmax(z))] == "C":
            s = s[~top]
        z = s.positions[:, 2]
        if s.get_chemical_symbols()[int(np.argmin(z))] == "C":
            s = s[z > z.min() + 0.3]
        s.positions[:, 2] -= s.positions[:, 2].min()
        return s
    if name in ("GaN", "AlN"):
        from phonons import relax
        b = relax(CRYSTALS[name]["build"]())
        s = surface(b, (0, 0, 1), nlayers, vacuum=0.0)
        P = np.array(size) if np.ndim(size) == 2 else np.diag([size[0], size[1], 1])
        return make_supercell(s, P)
    raise ValueError(name)


def bulk_reference(name, normal=None):
    from model import Side
    s = Side(name, normal)
    d = s.n_atoms ** (-1 / 3)
    return dict(K_n=s.rho * s.vL ** 2 / d, K_t=s.rho * s.vT ** 2 / d)


def Kref(a, b, na=None, nb=None):
    ra, rb = bulk_reference(a, na), bulk_reference(b, nb)
    hm = lambda x, y: 2 * x * y / (x + y)
    return dict(K_n=hm(ra["K_n"], rb["K_n"]), K_t=hm(ra["K_t"], rb["K_t"]))


def load_cache():
    if os.path.exists(CACHE):
        return json.load(open(CACHE))
    return {}


def save(key, val):
    c = load_cache(); c[key] = val
    json.dump(c, open(CACHE, "w"), indent=1)


JOBS = {}


def job(f):
    JOBS[f.__name__] = f
    return f


@job
def Au_PMMA(): return K_molecular_layer("Au")
@job
def Au_PMMA_rx(): return K_molecular_layer("Au", seeds=(0, 1), reactive=True)
@job
def Ti_PMMA(): return K_molecular_layer("Ti")
@job
def Ti_PMMA_rx(): return K_molecular_layer("Ti", seeds=(0, 1), reactive=True)
@job
def Ni_PMMA(): return K_molecular_layer("Ni")
@job
def Ni_PMMA_rx(): return K_molecular_layer("Ni", seeds=(0, 1), reactive=True)
@job
def PMMA_SiO2(): return K_molecular_layer("SiO2")
@job
def PMMA_SiO2_rx(): return K_molecular_layer("SiO2", seeds=(0, 1), reactive=True)


@job
def Al_Al():
    """Validation of the rigid-shift protocol against K_bulk = rho v^2/d:
    a bulk Al(111) crystal cut into two halves (registry kept)."""
    from ase.build import fcc111
    s = fcc111("Al", size=(2, 2, 12), a=fcc_a("Al"), vacuum=0.0, orthogonal=False)
    s.pbc = (True, True, False)
    z = s.positions[:, 2]
    up = z > np.median(z)
    A = np.linalg.norm(np.cross(s.cell[0], s.cell[1]))
    s.cell[2] = [0, 0, z.max() + 20]
    r = _scan(s, up, A)
    r["ref_rho_v2_over_d"] = bulk_reference("Al", "111")
    return r


@job
def Al_Si():
    lo = crystal_slab("Si", (3, 3), 6)          # 3 bilayers
    up = crystal_slab("Al", (4, 4), 4)
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 2.5); r["strain_upper"] = strain; return r


@job
def Al_GaN():
    P = np.array([[2, 1, 0], [-1, 1, 0], [0, 0, 1]])  # sqrt3 x sqrt3 R30
    lo = crystal_slab("GaN", P, 2)
    up = crystal_slab("Al", (2, 2), 4)
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 2.3); r["strain_upper"] = strain; return r


@job
def Cu_C_bare():
    lo = crystal_slab("C", (2, 2), 6)
    up = crystal_slab("Cu", (2, 2), 4)
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 2.0); r["strain_upper"] = strain; return r


@job
def Cu_C_H():
    from ase import Atoms
    lo = crystal_slab("C", (2, 2), 6)
    z = lo.positions[:, 2]
    top = lo[z > z.max() - 0.3]
    H = Atoms(["H"] * len(top), positions=top.positions + [0, 0, 1.10])
    lo = lo + H
    up = crystal_slab("Cu", (2, 2), 4)
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 2.2); r["strain_upper"] = strain; return r


@job
def Cu_TiC():
    lo = crystal_slab("TiC", (5, 5), 3)   # 5 x 3.055 = 15.27 A
    up = crystal_slab("Cu", (6, 6), 4)    # 6 x 2.569 = 15.41 A (-0.9%)
    s = lo.positions[:, 2]
    top_sym = lo.get_chemical_symbols()[int(np.argmax(s))]
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 2.1); r["strain_upper"] = strain; r["TiC_termination"] = top_sym; return r


@job
def TiC_C():
    lo = crystal_slab("C", (6, 6), 6)     # 6 x 2.517 = 15.10 A
    up = crystal_slab("TiC", (5, 5), 3)   # 5 x 3.055 = 15.27 A (-1.1%)
    up, strain = match(up, lo.cell)
    r = K_coherent(lo, up, 1.9); r["strain_upper"] = strain; return r


if __name__ == "__main__":
    names = sys.argv[1:] or list(JOBS)
    for nm in names:
        print(f"== {nm}", flush=True)
        r = JOBS[nm]()
        save(nm, r)
        print(json.dumps({k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in r.items()
                          if not isinstance(v, (list, dict))}), flush=True)
