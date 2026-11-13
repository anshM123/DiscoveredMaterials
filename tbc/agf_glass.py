"""Atomistic Green's function for crystal / amorphous-SiO2 interfaces with
MACE-OMAT-0 force constants (reference for the fast tier's glass treatment).

Superlattice [crystal x 4 periods | a-SiO2 x 3 identical copies] in the
11.45 A hexagonal cell of glass.py; positions and c relaxed with MACE (float32
FIRE + float64 polish); force constants with float32 MACE (tested: FC error
~1e-4 relative). Principal layers: one stacking period of the crystal, one
glass copy. Left lead = crystal, right lead = periodic glass.
Results cached in data/agf/<crystal>_aSiO2.json.
"""
import os, sys, json, warnings
import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import agf as A

SLABS = {"Al": ("fcc111", (4, 4, 3)), "Si": ("diamond111", (3, 3, 6))}


def crystal_period(name):
    from ase.build import fcc111, diamond111
    from bonding import fcc_a
    kind, size = SLABS[name]
    f = fcc111 if kind == "fcc111" else diamond111
    return f(name, size=size, a=fcc_a(name), orthogonal=False, periodic=True)


def build(name, glass, n_cryst=4, n_glass=3, gap=2.0):
    from ase import Atoms
    from bonding import match
    u = crystal_period(name)
    u, strain = match(u, glass.cell[:2])           # crystal strained onto the glass cell
    u.positions[:, 2] -= u.positions[:, 2].min()
    cz = u.cell[2, 2]; gz = glass.cell[2, 2]
    syms, pos, tag = [], [], []
    for k in range(n_cryst):
        for s, p in zip(u.get_chemical_symbols(), u.positions):
            syms.append(s); pos.append(p + [0, 0, k * cz]); tag.append(("L", k))
    z0 = n_cryst * cz - (cz - u.positions[:, 2].max()) + gap
    g = glass.copy(); g.wrap(); gp = g.positions.copy(); gp[:, 2] -= gp[:, 2].min()
    for k in range(n_glass):
        for s, p in zip(g.get_chemical_symbols(), gp):
            syms.append(s); pos.append(p + [0, 0, z0 + k * gz]); tag.append(("R", k))
    ctot = z0 + n_glass * gz + gap - (gz - gp[:, 2].max())
    cell = np.array([glass.cell[0], glass.cell[1], [0, 0, ctot]])
    return Atoms(syms, positions=pos, cell=cell, pbc=True), tag, strain


def _calcs():
    """MACE-OMAT-0 without D3 (covalent oxide/semiconductor/metal bonding;
    D3 makes MD/FC ~3x slower for no meaningful change here)."""
    from mace.calculators import mace_mp
    return (lambda: mace_mp(model="medium-omat-0", default_dtype="float32", device="cpu"),
            lambda: mace_mp(model="medium-omat-0", default_dtype="float64", device="cpu"))


def relax(sl):
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    calc32, calc = _calcs()
    s = sl.copy()
    s.calc = calc32(); FIRE(FrechetCellFilter(s, mask=[0, 0, 1, 0, 0, 0]), logfile=None).run(fmax=0.03, steps=400)
    s.calc = calc(); FIRE(FrechetCellFilter(s, mask=[0, 0, 1, 0, 0, 0]), logfile=None).run(fmax=0.02, steps=100)
    return s


def force_constants32(sl):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    calc32 = _calcs()[0]
    u = PhonopyAtoms(symbols=sl.get_chemical_symbols(), cell=sl.cell[:], scaled_positions=sl.get_scaled_positions())
    ph = Phonopy(u, supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")
    # one-sided displacements (no symmetry: 3N instead of 6N force calls)
    ph.generate_displacements(distance=0.01, is_plusminus=os.environ.get("AGF_PLUSMINUS", "0") == "1")
    c = calc32(); F = []
    n = len(ph.supercells_with_displacements)
    for k, sc in enumerate(ph.supercells_with_displacements):
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True); a.calc = c
        f = a.get_forces(); F.append(f - f.mean(axis=0))
        if k % 100 == 0:
            print(f"    displacement {k + 1}/{n}", flush=True)
    ph.forces = F; ph.produce_force_constants(); ph.symmetrize_force_constants()
    return ph


class GlassAGF(A.AGF):
    """L side: nL crystal periods; R side: nR glass copies (nR = 3).
    Device: L PL(nL-2), L PL(nL-1), R PL0. Right lead: R PL1 (H00) -> R PL2 (H01)."""

    def __init__(self, ph, tag, nL, nR):
        super().__init__(ph, tag, nL, (1, 1))
        self.nL, self.nR = nL, nR

    def lead_and_device(self, k):
        L = [(0, p) for p in range(self.nL)]
        R = [(1, p) for p in range(self.nR)]
        H00L, H01L = self.get(L[1], L[1], k), self.get(L[1], L[2], k)
        H00R, H01R = self.get(R[1], R[1], k), self.get(R[1], R[2], k)
        dev = [L[self.nL - 2], L[self.nL - 1], R[0]]
        sizes = [3 * len(self.idx[d]) for d in dev]
        N = sum(sizes); off = np.cumsum([0] + sizes)
        Hd = np.zeros((N, N), complex)
        for a_, da in enumerate(dev):
            for b_, db in enumerate(dev):
                Hd[off[a_]:off[a_ + 1], off[b_]:off[b_ + 1]] = self.get(da, db, k)
        VL = self.get(L[self.nL - 3], L[self.nL - 2], k)
        VR = self.get(R[0], R[1], k)
        return H00L, H01L, H00R, H01R, Hd, VL, VR, sizes


def run(name, nk=4, nw=200, glass_file=None):
    from ase.io import read
    from phonopy import load as phload
    os.makedirs(A.CACHE, exist_ok=True)
    glass = read(glass_file or os.path.join(HERE, "data", "glass", "aSiO2_32.extxyz"))
    label = f"{name}111_aSiO2"
    sl, tag, strain = build(name, glass)
    print(f"[{label}] SL atoms {len(sl)}, strain of crystal {strain:+.3f}", flush=True)
    fcfile = os.path.join(A.CACHE, f"{label}_ph.yaml")
    if os.path.exists(fcfile):
        ph = phload(fcfile, log_level=0)
    else:
        sl = relax(sl)
        print("    relaxed", flush=True)
        ph = force_constants32(sl)
        ph.save(fcfile, settings={"force_constants": True})
    G = GlassAGF(ph, tag, 4, 3)
    wmax = 2 * np.pi * 40e12 / A.OMEGA_UNIT
    wgrid = np.linspace(wmax / nw / 2, wmax, nw)
    K = A.kmesh(G.rec2d, nk)
    Xi = np.zeros(nw)
    for ik, k in enumerate(K):
        Xi += G.transmission(wgrid, k)
        print(f"    k {ik + 1}/{len(K)}", flush=True)
    Xi /= len(K)
    Ts = np.array([50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800], float)
    Gt = A.conductance(wgrid, Xi, G.A, Ts)
    res = dict(label=label, strain=strain, nk=nk, f_THz=list(wgrid * A.OMEGA_UNIT / 2 / np.pi / 1e12),
               Xi=list(Xi), A=G.A, T=list(Ts), G=list(Gt))
    json.dump(res, open(os.path.join(A.CACHE, f"{label}.json"), "w"), indent=1)
    print(json.dumps({"label": label, "G(80K)": round(float(Gt[1]) / 1e6, 1), "G(300K)": round(float(Gt[6]) / 1e6, 1)}), flush=True)
    return res


if __name__ == "__main__":
    run(sys.argv[1])
