"""Atomistic Green's function (AGF) phonon transmission across a coherent
interface, with MACE-OMAT-0 force constants.

This is the slow, rigorous reference ("atomistic ground truth" for a perfect,
harmonic, lattice-matched interface) against which the fast model in model.py
is checked. It is NOT part of the fast predictor.

Method
------
1. Oriented periodic units of both crystals along the interface normal
   (ase.build.surface, periodic=True). The R unit is strained in-plane onto
   the L in-plane lattice (lattice-matched pairs only, misfit <= ~4%).
2. Superlattice [L x nL | R x nR], periodic; positions and the c-length are
   relaxed with MACE-OMAT-0.
3. Harmonic force constants with phonopy on an (n x n x 1) in-plane supercell.
4. Atoms are grouped into principal layers (PL) of >= PL_MIN Angstrom along
   the normal. Leads: bulk-like PLs inside each material; device: the 2 PLs
   on each side of one interface. k_par-resolved blocks are built from the
   supercell force constants with minimal-image in-plane phases.
5. Sancho-Rubio decimation for lead surface Green's functions;
   Xi(w,k) = Tr[Gamma_L G Gamma_R G^+];
   G = 1/(A N_k) sum_k int dw/2pi  hbar w Xi(w,k) dn/dT.
Same Landauer convention (emitted-phonon temperatures) as model.py.
"""
import os, sys, json, warnings
import numpy as np

# silence third-party library noise only; this package's own warnings
# (missing e-ph data, estimated g, fallbacks) must reach the user
for _cat in (DeprecationWarning, FutureWarning):
    warnings.filterwarnings("ignore", category=_cat)
warnings.filterwarnings("ignore", module=r"(torch|e3nn|mace|ase|phonopy|spglib|numpy|torch_dftd)(\.|$)")
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "data", "agf")
HBAR = 1.054571817e-34
KB = 1.380649e-23
OMEGA_UNIT = 9.82269385e13       # sqrt(eV / (Angstrom^2 amu)) in rad/s
PL_MIN = 8.5                     # Angstrom, principal-layer thickness


def oriented_unit(name, miller):
    """Relaxed bulk of `name`, re-cut so that a3 spans one period along the
    normal and a1, a2 are in-plane. Miller indices refer to the spglib
    standardised conventional cell (cubic for fcc/bcc/diamond/zincblende/
    rocksalt, hexagonal for wurtzite/hcp/corundum ...), so any crystal works."""
    from ase.build import surface
    from ase import Atoms
    import spglib
    from phonons import relax
    from materials import CRYSTALS
    b = relax(CRYSTALS[name]["build"]())
    lat, frac, nums = spglib.standardize_cell((b.cell[:], b.get_scaled_positions(), b.numbers),
                                              to_primitive=False, no_idealize=False, symprec=1e-3)
    conv = Atoms(numbers=nums, cell=lat, scaled_positions=frac, pbc=True)
    digits = [int(c) for c in miller.replace("-", "")]
    if len(digits) == 4:                 # hkil -> hkl
        digits = [digits[0], digits[1], digits[3]]
    if tuple(digits) == (0, 0, 1) and abs(conv.cell.angles()[2] - 120) < 1:
        u = conv.copy()                  # hexagonal c-axis: the conventional cell already is the unit
    else:
        u = surface(conv, tuple(digits), 1, vacuum=None, periodic=True)
    # rotate so that a1 || x and the in-plane normal || z
    c = u.cell.array
    n = np.cross(c[0], c[1]); n /= np.linalg.norm(n)
    e1 = c[0] / np.linalg.norm(c[0]); e2 = np.cross(n, e1)
    R = np.array([e1, e2, n])
    u.set_cell(c @ R.T, scale_atoms=True)
    if u.cell[2, 2] < 0:
        raise RuntimeError("left-handed oriented cell")
    return u


def superlattice(L, nL, R, nR):
    """Stack nL periods of L and nR periods of R (R strained in-plane to L)."""
    from ase import Atoms
    cL, cR = L.cell.array, R.cell.array
    # in-plane strain of R onto L: map R's (a1,a2) to L's (a1,a2)
    F = np.linalg.solve(cR[:2, :2], cL[:2, :2])      # 2x2 in-plane map (row-vector convention)
    strain = float(np.sqrt(abs(np.linalg.det(F))) - 1)
    def map_R(v):
        out = v.copy()
        out[..., :2] = v[..., :2] @ F
        return out
    cR3 = map_R(cR[2])
    syms, pos, tag = [], [], []
    for k in range(nL):
        for s, p in zip(L.get_chemical_symbols(), L.positions):
            syms.append(s); pos.append(p + k * cL[2]); tag.append(("L", k))
    off = nL * cL[2]
    Rpos = map_R(R.positions)
    for k in range(nR):
        for s, p in zip(R.get_chemical_symbols(), Rpos):
            syms.append(s); pos.append(off + p + k * cR3); tag.append(("R", k))
    cell = np.array([cL[0], cL[1], off + nR * cR3])
    sl = Atoms(syms, positions=pos, cell=cell, pbc=True)
    return sl, tag, strain


def relax_sl(sl):
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS
    from phonons import calculator
    s = sl.copy(); s.calc = calculator()
    BFGS(FrechetCellFilter(s, mask=[0, 0, 1, 0, 0, 0]), logfile=None).run(fmax=2e-3, steps=400)
    return s


def force_constants(sl, nsc):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    from phonons import calculator
    u = PhonopyAtoms(symbols=sl.get_chemical_symbols(), cell=sl.cell[:], scaled_positions=sl.get_scaled_positions())
    ph = Phonopy(u, supercell_matrix=np.diag([nsc, nsc, 1]), primitive_matrix="P")
    ph.generate_displacements(distance=0.01)
    calc = calculator()
    F = []
    for k, sc in enumerate(ph.supercells_with_displacements):
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        a.calc = calc
        f = a.get_forces(); F.append(f - f.mean(axis=0))
        if k % 20 == 0:
            print(f"    displacement {k+1}/{len(ph.supercells_with_displacements)}", flush=True)
    ph.forces = F
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    return ph


class AGF:
    def __init__(self, ph, tag, nPL_side, periods_per_PL, fc_tol=1e-6):
        sc = ph.supercell
        self.masses = np.array(ph.primitive.masses)
        fc = ph.force_constants
        p2s = np.asarray(ph.primitive.p2s_map)
        s2p = np.asarray(ph.primitive.s2p_map)
        npr = len(p2s)
        fc_rows = fc if fc.shape[0] == npr else fc[p2s]
        s_to_prim = {int(s): i for i, s in enumerate(p2s)}
        prim_of_sc = np.array([s_to_prim[int(t)] for t in s2p])
        pos = sc.positions; latt = sc.cell; inv = np.linalg.inv(latt)
        self.z = ph.primitive.positions[:, 2]
        self.cz = ph.primitive.cell[2, 2]
        ppl = periods_per_PL if isinstance(periods_per_PL, (tuple, list)) else (periods_per_PL, periods_per_PL)
        self.pl = [(0 if s_ == "L" else 1, k // ppl[0 if s_ == "L" else 1]) for s_, k in tag]
        self.nPL = nPL_side
        self.idx = {}
        for a_, key in enumerate(self.pl):
            self.idx.setdefault(key, []).append(a_)
        self.idx = {k: np.array(v) for k, v in self.idx.items()}
        loc = {}
        for key, at in self.idx.items():
            for r, a_ in enumerate(at):
                loc[int(a_)] = (key, r)
        # group all non-negligible force-constant pairs by (PL_i, PL_j, z-image)
        groups = {}
        # all periodic images within one supercell in each direction; a force
        # constant is shared equally among equidistant minimal images
        # (multiplicity, as in phonopy) so that D(k) is exactly Hermitian even
        # for even supercells, where pairs lie exactly half a cell apart
        shifts = np.array([[a_, b_, c_] for a_ in (-1, 0, 1) for b_ in (-1, 0, 1) for c_ in (-1, 0, 1)], float)
        for i in range(npr):
            d0 = pos - pos[p2s[i]]
            fr = d0 @ inv; fr -= np.round(fr)
            norms = np.abs(fc_rows[i]).max(axis=(1, 2))
            for sj in np.nonzero(norms > fc_tol)[0]:
                cand = (fr[sj] + shifts) @ latt
                dist = np.linalg.norm(cand, axis=1)
                sel = cand[dist < dist.min() + 1e-4]
                j = int(prim_of_sc[sj])
                ki, ri = loc[i]; kj, rj = loc[j]
                f = fc_rows[i, sj] / np.sqrt(self.masses[i] * self.masses[j]) / len(sel)
                for d in sel:
                    img = int(np.round((self.z[i] + d[2] - self.z[j]) / self.cz))
                    g = groups.setdefault((ki, kj, img), ([], [], [], []))
                    g[0].append(ri); g[1].append(rj); g[2].append(d[:2]); g[3].append(f)
        self.groups = {k: (np.array(v[0]), np.array(v[1]), np.array(v[2]), np.array(v[3])) for k, v in groups.items()}
        self.A = abs(np.linalg.det(ph.primitive.cell[:2, :2]))     # Angstrom^2
        self.rec2d = 2 * np.pi * np.linalg.inv(ph.primitive.cell[:2, :2]).T

    def get(self, a, b, k):
        """Mass-weighted block between principal layers a and b at k_par."""
        na, nb = 3 * len(self.idx[a]), 3 * len(self.idx[b])
        B = np.zeros((na // 3, nb // 3, 3, 3), complex)
        g = self.groups.get((a, b, 0))
        if g is not None:
            ri, rj, d, f = g
            ph = np.exp(1j * (d @ k))
            np.add.at(B, (ri, rj), f * ph[:, None, None])
        return B.transpose(0, 2, 1, 3).reshape(na, nb)

    def lead_and_device(self, k):
        """Left lead: L PL1 (H00) -> L PL2 (H01). Device: L PL2, L PL3, R PL0, R PL1.
        Right lead: R PL2 (H00) -> R PL3 (H01)."""
        nP = self.nPL
        L = [(0, p) for p in range(nP)]
        Rr = [(1, p) for p in range(nP)]
        H00L, H01L = self.get(L[1], L[1], k), self.get(L[1], L[2], k)
        H00R, H01R = self.get(Rr[nP - 2], Rr[nP - 2], k), self.get(Rr[nP - 2], Rr[nP - 1], k)
        dev = [L[nP - 2], L[nP - 1], Rr[0], Rr[1]]
        sizes = [3 * len(self.idx[d]) for d in dev]
        N = sum(sizes); off = np.cumsum([0] + sizes)
        Hd = np.zeros((N, N), complex)
        for a_, da in enumerate(dev):
            for b_, db in enumerate(dev):
                Hd[off[a_]:off[a_ + 1], off[b_]:off[b_ + 1]] = self.get(da, db, k)
        VL = self.get(L[nP - 3], L[nP - 2], k)     # lead surface cell (L PL1 type) -> device PL0
        VR = self.get(Rr[1], Rr[2], k)             # device last PL -> lead surface cell (R PL2 type)
        return H00L, H01L, H00R, H01R, Hd, VL, VR, sizes

    @staticmethod
    def surface_gf(w2, H00, H01, left, eta=1e-6, tol=1e-10, itmax=200):
        """Sancho-Rubio. Returns surface GF of the semi-infinite lead whose
        surface faces the device (right end of a left lead, left end of a right lead)."""
        n = H00.shape[0]
        I = np.eye(n)
        E = (w2 + 1j * eta) * I
        if left:
            a, b = H01.conj().T, H01        # couplings towards the bulk: cell n -> n-1 = H10
        else:
            a, b = H01, H01.conj().T
        es = H00.copy(); e = H00.copy()
        for _ in range(itmax):
            g = np.linalg.inv(E - e)
            agb = a @ g @ b; bga = b @ g @ a
            es = es + agb
            e = e + agb + bga
            a = a @ g @ a; b = b @ g @ b
            if np.abs(a).max() < tol and np.abs(b).max() < tol:
                break
        return np.linalg.inv(E - es)

    def transmission(self, wgrid, k, eta_rel=1e-3):
        H00L, H01L, H00R, H01R, Hd, VL, VR, sizes = self.lead_and_device(k)
        N = Hd.shape[0]; nl = sizes[0]; nr = sizes[-1]
        Xi = np.zeros(len(wgrid))
        for t, w in enumerate(wgrid):
            w2 = w * w
            eta = eta_rel * max(w2, 1e-4) + 1e-8
            gL = self.surface_gf(w2, H00L, H01L, left=True, eta=eta)
            gR = self.surface_gf(w2, H00R, H01R, left=False, eta=eta)
            SL = VL.conj().T @ gL @ VL
            SR = VR @ gR @ VR.conj().T
            Sig = np.zeros((N, N), complex)
            Sig[:nl, :nl] = SL
            Sig[N - nr:, N - nr:] = SR
            Gd = np.linalg.inv((w2 + 1j * eta) * np.eye(N) - Hd - Sig)
            GamL = 1j * (SL - SL.conj().T)
            GamR = 1j * (SR - SR.conj().T)
            Xi[t] = np.real(np.trace(GamL @ Gd[:nl, N - nr:] @ GamR @ Gd[:nl, N - nr:].conj().T))
        return Xi


def kmesh(rec2d, n):
    f = (np.arange(n) + 0.5) / n - 0.5
    K = np.array([[a, b] for a in f for b in f])
    return K @ rec2d


def conductance(wgrid, Xi_avg, A_ang2, Ts):
    """wgrid in units of OMEGA_UNIT; Xi averaged over k; A in Angstrom^2."""
    w = wgrid * OMEGA_UNIT
    dw = w[1] - w[0]
    out = []
    for T in Ts:
        x = HBAR * w / (KB * T)
        dndT = np.where(x < 700, x / T * np.exp(np.minimum(x, 700)) / np.expm1(np.minimum(x, 700)) ** 2, 0.0)
        out.append(np.sum(HBAR * w * Xi_avg * dndT) * dw / (2 * np.pi) / (A_ang2 * 1e-20))
    return np.array(out)


def run(nameL, mL, nameR, mR, nsc=0, nk=8, nw=160, periods=None, label=None, inplane_min=12.0):
    os.makedirs(CACHE, exist_ok=True)
    label = label or f"{nameL}{mL}_{nameR}{mR}"
    uL, uR = oriented_unit(nameL, mL), oriented_unit(nameR, mR)
    ppl_L = max(1, int(np.ceil(PL_MIN / uL.cell[2, 2])))
    ppl_R = max(1, int(np.ceil(PL_MIN / uR.cell[2, 2])))
    nP = 4
    sl, tag, strain = superlattice(uL, nP * ppl_L, uR, nP * ppl_R)
    print(f"[{label}] SL atoms {len(sl)}, periods/PL L {ppl_L} R {ppl_R}, in-plane strain of R {strain:+.3f}", flush=True)
    if not nsc:
        nsc = int(np.ceil(inplane_min / min(np.linalg.norm(uL.cell[0]), np.linalg.norm(uL.cell[1]))))
    print(f"    in-plane FC supercell {nsc}x{nsc}: {len(sl) * nsc * nsc} atoms", flush=True)
    fcfile = os.path.join(CACHE, f"{label}_ph.yaml")
    if os.path.exists(fcfile):
        from phonopy import load as phload
        ph = phload(fcfile, log_level=0)
        print("    loaded cached force constants", flush=True)
    else:
        sl = relax_sl(sl)
        ph = force_constants(sl, nsc)
        ph.save(fcfile, settings={"force_constants": True})
    agf = AGF(ph, tag, nP, (ppl_L, ppl_R))
    # frequency grid up to max phonon frequency of both bulk materials
    from phonons import load
    fmax = max(load(nameL)["f_max"], load(nameR)["f_max"]) * 1.03
    wmax = 2 * np.pi * fmax * 1e12 / OMEGA_UNIT
    wgrid = np.linspace(wmax / nw / 2, wmax, nw)
    K = kmesh(agf.rec2d, nk)
    Xi = np.zeros(nw)
    for ik, k in enumerate(K):
        Xi += agf.transmission(wgrid, k)
        if ik % 8 == 0:
            print(f"    k {ik+1}/{len(K)}", flush=True)
    Xi /= len(K)
    Ts = np.array([50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800], float)
    G = conductance(wgrid, Xi, agf.A, Ts)
    res = dict(label=label, L=nameL, R=nameR, normal_L=mL, normal_R=mR, strain_R=strain,
               nsc=nsc, nk=nk, f_THz=list(wgrid * OMEGA_UNIT / 2 / np.pi / 1e12), Xi=list(Xi),
               A=agf.A, T=list(Ts), G=list(G))
    json.dump(res, open(os.path.join(CACHE, f"{label}.json"), "w"), indent=1)
    print(json.dumps({"label": label, "G(300K) MW/m2K": round(float(np.interp(300, Ts, G)) / 1e6, 1)}), flush=True)
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    kw = {}
    for x in a[4:]:
        k_, v = x.split("=")
        kw[k_] = v if k_ == "label" else (float(v) if "." in v else int(v))
    run(a[0], a[1], a[2], a[3], **kw)
