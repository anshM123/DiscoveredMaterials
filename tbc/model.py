"""Thermal boundary conductance: baselines (Debye AMM / DMM) and the
transferable model described in README.md.

All conductances are returned in W m^-2 K^-1 unless stated otherwise.
Angular frequency grid W (rad/s) is shared with phonons.py.
"""
import os
import numpy as np
from phonons import W, DW, load, T_GRID
from materials import CRYSTALS, AMORPHOUS, METALS

HBAR = 1.054571817e-34
KB = 1.380649e-23
AMU = 1.66053906660e-27
THETA_SIDE = "stiff"
INELASTIC_AMORPHOUS = False   # inelastic channel only for crystal/crystal interfaces
K_RESOLVED_SPECULAR = True    # match specular channels at equal k_par
SPECTRAL_IMPEDANCE = True     # bonding factor with lattice impedance rho*v_g(w)
SMOOTH_THZ = 0.0             # Gaussian smoothing (THz) of histogrammed crystal spectra before nonlinear matching
LOWW_DEBYE = False           # (H12, tested, not adopted) exact long-wavelength (Debye) M, D below LOWW_DEBYE_MESH q-mesh spacings
LOWW_DEBYE_MESH = 3.0
LOWW_CONTINUUM = False       # (H12, tested, not adopted) tau_spec held at its value at that frequency below it (continuum limit)
K_IMPEDANCE = True            # impedance factor per (k_par, w) in the specular channel
SPEC_EFFICIENCY = 1.0         # specular efficiency between different crystals (set from params.json)
SPEC_EFFICIENCY_HETERO = None # same for heterovalent III-V/IV pairs (None: use SPEC_EFFICIENCY)
W_DIFFUSON = 1.0              # interfacial coupling weight of glass diffusons (set from params.json)
_GLASS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "glass")
GLASS_FLUX = {"SiO2": os.path.join(_GLASS_DIR, "aSiO2_32_flux.npz"),       # glass.py SiO2
              "aAl2O3": os.path.join(_GLASS_DIR, "aAl2O3_20_flux.npz")}   # glass.py aAl2O3
try:                          # frozen global parameters
    import json as _json, os as _os
    _P = _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "params.json")))
    SPEC_EFFICIENCY = float(_P.get("spec_efficiency", SPEC_EFFICIENCY))
    if "spec_efficiency_hetero" in _P:
        SPEC_EFFICIENCY_HETERO = float(_P["spec_efficiency_hetero"])
    W_DIFFUSON = float(_P.get("w_diffuson", W_DIFFUSON))
except Exception:             # pragma: no cover
    pass
NK_SPEC = 16                  # k_par grid for channel_spectra
NZ_SPEC = 32                  # q_z points per column (flux captured ~95%, G converged to ~1%)


def dfdT(w, T):
    """d n_BE / dT."""
    x = np.clip(HBAR * w / (KB * T), 1e-12, 700)
    ex = np.exp(x)
    return x / T * ex / (ex - 1.0) ** 2


def bose(w, T):
    x = np.clip(HBAR * w / (KB * T), 1e-12, 700)
    return 1.0 / np.expm1(x)


class Side:
    """One side of an interface: spectral flux M(w) normal to the interface."""

    def __init__(self, name, normal=None):
        self.name = name
        if name in CRYSTALS:
            d = load(name)
            if normal is None:
                normal = next(iter(CRYSTALS[name]["normals"]))
            self.normal = normal
            self.M = d[f"M_{normal}"]
            self.Mperp = d[f"Mperp_{normal}"]
            self.D = d["D"]
            # 2D Bravais lattice of the interface plane (Angstrom), any normal
            nvec = np.asarray(CRYSTALS[name]["normals"][normal], float)
            self.plane = plane_lattice(d["a_cell"], nvec)
            self.a_surf = float(np.linalg.norm(self.plane[0]))
            # MACE (PBE-like) lattice constants carry errors of 1-2%, which are
            # comparable to the misfits that decide coherence. The coherence
            # test therefore uses the plane rescaled to the input (experimental)
            # atomic volume: lat_scale = (V_input / V_MACE)^(1/3) per atom.
            try:
                ref = CRYSTALS[name]["build"]()
                self.lat_scale = float((ref.get_volume() / len(ref) * 1e-30 / (d["V"] / d["nat"])) ** (1 / 3))
            except Exception:
                self.lat_scale = 1.0
            self.rho = float(d["rho"])
            self.vL, self.vT = float(d["vL"]), float(d["vT"])
            self.n_atoms = float(d["n_atoms"])
            self.m_avg = float(d["m_avg"])
            self.crystal = True
            # Long-wavelength limit. The q-mesh resolves wavevectors only down to
            # ~|b_i|/N_i, so below a few mesh spacings the histogrammed M(w), D(w)
            # are spiky (e.g. AlN: M = 0 at 0.5 THz). There the exact continuum
            # (Debye) forms hold: D = w^2/(2 pi^2) (1/vL^3 + 2/vT^3),
            # M = w^2/(8 pi^2) (1/vL^2 + 2/vT^2), with the MACE sound velocities.
            if SMOOTH_THZ > 0:
                self.M, self.Mperp, self.D = (_smooth(x, SMOOTH_THZ) for x in (self.M, self.Mperp, self.D))
            self.w_debye = 0.0
            if LOWW_DEBYE:
                B = 2 * np.pi * np.linalg.inv(np.asarray(d["a_cell"], float) * 1e-10).T
                mesh = CRYSTALS[name].get("mesh", [24, 24, 24])
                qres = min(np.linalg.norm(B[i]) / mesh[i] for i in range(3))
                self.w_debye = LOWW_DEBYE_MESH * self.vT * qres
                low = W < self.w_debye
                ML = W ** 2 / (8 * np.pi ** 2 * self.vL ** 2); MT = 2 * W ** 2 / (8 * np.pi ** 2 * self.vT ** 2)
                self.M = np.where(low, ML + MT, self.M)
                self.Mperp = np.where(low, 0.5 * ML + 0.5 * MT, self.Mperp)
                self.D = np.where(low, W ** 2 / (2 * np.pi ** 2) * (1 / self.vL ** 3 + 2 / self.vT ** 3), self.D)
        else:
            a = AMORPHOUS[name]
            self.normal = None
            self.rho, self.vL, self.vT = a.rho, a.vL, a.vT
            self.n_atoms, self.m_avg = a.n_atoms, a.m_avg
            self.crystal = False
            ML = np.zeros_like(W); MT = np.zeros_like(W)
            for v, g, arr in [(self.vL, 1, ML), (self.vT, 2, MT)]:
                wc = v * self.kD
                D = g * W ** 2 / (2 * np.pi ** 2 * v ** 3)
                arr[:] = np.where(W <= wc, D * v / 4, 0.0)
            self.M = ML + MT
            self.D = np.where(W <= self.vL * self.kD, W ** 2 / (2 * np.pi ** 2 * self.vL ** 3), 0.0) + \
                np.where(W <= self.vT * self.kD, 2 * W ** 2 / (2 * np.pi ** 2 * self.vT ** 3), 0.0)
            self.a_surf = None
            self.plane = None
            # isotropic solid: flux-weighted fraction of displacement along n is 1/2
            self.Mperp = 0.5 * ML + 0.5 * MT
            # Glasses with an atomistic (MACE melt-quench) model: only a fraction of
            # the Debye channels are propagons; the rest are diffusons, which couple
            # across an interface with the global weight W_DIFFUSON:
            #   M_eff = M_prop + W_DIFFUSON * (M_Debye - M_prop),
            # M_prop = propagating flux of the periodic MACE glass (glass.py).
            gf = GLASS_FLUX.get(name)
            if gf is not None and os.path.exists(gf):
                Mp = np.load(gf)["M"]
                frac = np.where(self.M > 0, np.clip(np.divide(Mp, self.M, out=np.zeros_like(Mp), where=self.M > 0), 0, 1), 0)
                f = frac + W_DIFFUSON * (1 - frac)
                self.M = self.M * f
                self.Mperp = self.Mperp * f
        nz = np.nonzero(self.M > 1e-6 * self.M.max())[0]
        self.w_max = W[nz[-1]] + DW / 2

    @property
    def kD(self):
        return (6 * np.pi ** 2 * self.n_atoms) ** (1 / 3)

    def debye_branches(self):
        """(velocity, cutoff, multiplicity) of Debye branches (Swartz-Pohl)."""
        return [(self.vL, self.vL * self.kD, 1), (self.vT, self.vT * self.kD, 2)]

    def G_rad(self, T):
        """Phonon radiation limit (all phonons transmitted)."""
        return np.sum(HBAR * W * dfdT(W, T) * self.M) * DW

    def __repr__(self):
        return f"Side({self.name},{self.normal})"


# ---------------------------------------------------------------- baselines
def _debye_M(side):
    M = np.zeros_like(W)
    for v, wc, g in side.debye_branches():
        M += np.where(W <= wc, g * W ** 2 / (2 * np.pi ** 2 * v ** 3) * v / 4, 0.0)
    return M


def G_DMM_debye(s1, s2, T):
    """Diffuse mismatch model, Debye dispersion (Swartz & Pohl 1989)."""
    M1, M2 = _debye_M(s1), _debye_M(s2)
    tau = np.divide(M2, M1 + M2, out=np.zeros_like(W), where=(M1 + M2) > 0)
    return np.sum(HBAR * W * dfdT(W, T) * M1 * tau) * DW


def _amm_angle_avg(Z1, Z2, v1, v2, n=2000):
    """2 * int_0^{pi/2} tau(theta) cos(theta) sin(theta) dtheta, no mode conversion."""
    th = (np.arange(n) + 0.5) * (np.pi / 2) / n
    s2 = v2 / v1 * np.sin(th)
    ok = s2 < 1
    c1 = np.cos(th)
    c2 = np.sqrt(np.clip(1 - s2 ** 2, 0, 1))
    tau = np.where(ok, 4 * Z1 * Z2 * c1 * c2 / (Z1 * c2 + Z2 * c1) ** 2, 0.0)   # acoustic Fresnel (energy)
    return 2 * np.sum(tau * c1 * np.sin(th)) * (np.pi / 2) / n


def G_AMM_debye(s1, s2, T):
    """Acoustic mismatch model, Debye dispersion (Little 1959), evaluated from side 1.
    Each branch transmits only below both sides' Debye cutoffs (elastic
    scattering), which makes it reciprocal: G(1,2) = G(2,1)."""
    G = 0.0
    for (v1, wc1, g), (v2, wc2, _) in zip(s1.debye_branches(), s2.debye_branches()):
        if not (v1 > 0 and v2 > 0):
            continue
        a = _amm_angle_avg(s1.rho * v1, s2.rho * v2, v1, v2)
        M1 = np.where(W <= min(wc1, wc2), g * W ** 2 / (2 * np.pi ** 2 * v1 ** 3) * v1 / 4, 0.0)
        G += np.sum(HBAR * W * dfdT(W, T) * M1) * DW * a
    return G


# ---------------------------------------------------------------- the model
COHERENCE_STRAIN = 0.05   # max misfit for a (locally) coherent match: at 3-5% misfit, epitaxial
                          # interfaces relax by misfit dislocations ~a/f ~ 6-10 nm apart and are
                          # coherent in between (patches >> dominant phonon wavelengths)


def _smooth(y, sigma_thz, axis=-1):
    """Gaussian smoothing along the frequency grid (preserves the integral)."""
    if sigma_thz <= 0:
        return y
    sig = sigma_thz * 2 * np.pi * 1e12 / DW
    n = int(np.ceil(4 * sig))
    k = np.exp(-0.5 * (np.arange(-n, n + 1) / sig) ** 2); k /= k.sum()
    return np.apply_along_axis(lambda v: np.convolve(v, k, mode="same"), axis, y)


def spec_efficiency(a, b):
    """Specular efficiency for a pair of different crystals: channel counting cannot
    see eigenvector (polarisation/sublattice) mismatch. Heterovalent III-V/group-IV
    interfaces, whose electron-deficient interfacial bonds relax to soft, stretched
    ones (MACE-AGF structures, EXPERIMENTS), get their own value. Both are global
    constants calibrated on the MACE-AGF references (calibrate_spec.py)."""
    from materials import heterovalent
    if SPEC_EFFICIENCY_HETERO is not None and heterovalent(a, b) and "III-V" in heterovalent(a, b):
        return SPEC_EFFICIENCY_HETERO
    return SPEC_EFFICIENCY


def plane_lattice(cell, normal, nmax=4):
    """Reduced 2D basis (2x2, in-plane orthonormal coordinates, Angstrom) of
    the Bravais-lattice plane perpendicular to `normal`. Works for any
    orientation: (111), (100), (110), (0001), ... Frame: e1 along the shortest
    in-plane lattice vector, second vector at a positive angle."""
    cell = np.asarray(cell, float)
    n = normal / np.linalg.norm(normal)
    r = np.arange(-nmax, nmax + 1)
    I = np.array(np.meshgrid(r, r, r)).reshape(3, -1).T
    V = I @ cell
    L = np.linalg.norm(V, axis=1)
    ok = (L > 1e-6) & (np.abs(V @ n) < 1e-4 * np.maximum(L, 1))
    V, L = V[ok], L[ok]
    order = np.argsort(L, kind="stable")
    v1 = V[order[0]]
    for k in order[1:]:
        if np.linalg.norm(np.cross(v1, V[k])) > 1e-3 * L[order[0]] * L[k]:
            v2 = V[k]; break
    e1 = v1 / np.linalg.norm(v1); e2 = np.cross(n, e1)
    b = np.array([[v1 @ e1, v1 @ e2], [v2 @ e1, v2 @ e2]])
    if b[1, 1] < 0:
        b[1] = -b[1]
    return b


def _area(b):
    return abs(b[0, 0] * b[1, 1] - b[0, 1] * b[1, 0])


def coincidence_cell(s1, s2, nmax=12, tol=COHERENCE_STRAIN, dtheta=1.0):
    """Smallest common supercell of the two interface-plane lattices in the
    aligned orientation (shortest in-plane lattice vectors of both sides
    parallel: close-packed rows parallel for (111)/(0001), cube-on-cube for
    (100)/(110)), within misfit strain tol. Returns (Lc, p_inf) with Lc the
    longest common supercell vector (Angstrom) and p_inf = area of the smaller
    primitive plane cell / common cell area; None if a side is amorphous or no
    match exists (-> diffuse)."""
    if s1.plane is None or s2.plane is None:
        return None
    r = np.arange(-nmax, nmax + 1)
    I = np.array(np.meshgrid(r, r)).reshape(2, -1).T
    I = I[np.any(I != 0, axis=1)]

    def cand(b):
        U = I @ b
        return U, np.linalg.norm(U, axis=1), np.degrees(np.arctan2(U[:, 1], U[:, 0]))
    U, LU, TU = cand(s1.plane * getattr(s1, "lat_scale", 1.0))
    Wv, LW, TW = cand(s2.plane * getattr(s2, "lat_scale", 1.0))
    pairs = []
    for i in np.argsort(LU):
        if LU[i] > nmax * max(s1.a_surf, s2.a_surf) * 1.05:
            break
        m = (np.abs(np.log(LW / LU[i])) <= tol) & (np.abs((TW - TU[i] + 180) % 360 - 180) <= dtheta)
        for j in np.nonzero(m)[0]:
            pairs.append((LU[i], U[i], Wv[j]))
    best = None
    for a in range(len(pairs)):
        for c in range(a + 1, len(pairs)):
            u1, w1, u2, w2 = pairs[a][1], pairs[a][2], pairs[c][1], pairs[c][2]
            Au, Aw = abs(u1[0] * u2[1] - u1[1] * u2[0]), abs(w1[0] * w2[1] - w1[1] * w2[0])
            if min(Au, Aw) < 1e-6:
                continue
            A = max(Au, Aw)
            # the full 2x2 map between the two supercells must be a small
            # strain (singular values within tol) with negligible rotation;
            # rejects spurious cells spanned by long, nearly parallel vectors
            Tm = np.column_stack([w1, w2]) @ np.linalg.inv(np.column_stack([u1, u2]))
            Uu, sv, Vt = np.linalg.svd(Tm)
            Rm = Uu @ Vt
            if np.max(np.abs(np.log(sv))) > tol or abs(np.degrees(np.arctan2(Rm[1, 0], Rm[0, 0]))) > dtheta:
                continue
            Lc = max(np.linalg.norm(u1), np.linalg.norm(u2), np.linalg.norm(w1), np.linalg.norm(w2))
            if best is None or A < best[0] - 1e-6 or (abs(A - best[0]) < 1e-6 and Lc < best[1]):
                best = (A, Lc)
    if best is None:
        return None
    p_inf = min(_area(s1.plane) * getattr(s1, "lat_scale", 1.0) ** 2, _area(s2.plane) * getattr(s2, "lat_scale", 1.0) ** 2) / best[0]
    return best[1], min(1.0, p_inf)


def coincidence_length(s1, s2, **kw):
    c = coincidence_cell(s1, s2, **kw)
    return None if c is None else c[0]


def spring_factor(w, K, Z1, Z2):
    """Energy transmission of a plane wave through an interfacial spring K
    (N m^-3) relative to a welded contact."""
    zr = np.divide(Z1 * Z2, Z1 + Z2, out=np.zeros_like(np.asarray(Z1 + Z2, float)), where=np.asarray(Z1 + Z2) > 0)
    return 1.0 / (1.0 + (w * zr / K) ** 2)


class Interface:
    """Phonon-phonon interface between Side s1 and Side s2.

    Parameters
    ----------
    K : (K_normal, K_shear) interfacial stiffness per area of the actual
        interface (N m^-3) or None for an ideally (bulk-like) bonded interface.
    Kref : (K_normal, K_shear) stiffness of a bulk-like bond across the same
        plane; the bonding factor is spring(K)/spring(Kref), capped at 1.
    sharp : bool, whether the inelastic channel is open (atomically sharp
        interface). Intermixed/diffuse interfaces suppress it (Li et al. 2022).
    """

    def __init__(self, s1, s2, K=None, Kref=None, sharp=True, Lc="auto", epitaxial=True):
        """epitaxial: True for a coherent/epitaxial interface (specularity from
        the lattice coincidence); False for polycrystalline or sputtered films,
        which scatter diffusely (p = 0). A property of the specified interface
        structure, not a fitted parameter."""
        self.s1, self.s2 = s1, s2
        if not epitaxial:
            Lc = None
        self.K, self.Kref, self.sharp = K, Kref, sharp
        M1, M2 = s1.M, s2.M
        self.tau_dmm = np.divide(M2, M1 + M2, out=np.zeros_like(W), where=(M1 + M2) > 0)
        cc = coincidence_cell(s1, s2) if Lc == "auto" else None
        if Lc == "auto":
            self.Lc, self.p_inf = (cc if cc is not None else (None, 0.0))
        else:
            self.Lc = Lc
            self.p_inf = None if Lc is None else min(1.0, (min(s1.a_surf, s2.a_surf) / Lc) ** 2)
        self.tau_spec = self.specular_transmission()
        self.p_spec = self.specularity(W)
        self.tau_el = self.p_spec * self.tau_spec + (1 - self.p_spec) * self.tau_dmm
        self.S = self.bond_factor(W)

    def specular_transmission(self):
        """Momentum-conserving (specular) limit. Channels can only be matched at
        the same in-plane wavevector k_par (aligned orientation, k_par compared
        in the fractional coordinates of each side's interface-plane reciprocal
        cell), so the transmissible flux is sum_k min(M1_k, M2_k), times the
        impedance factor 4 Z1 Z2/(Z1+Z2)^2 (Z_i = rho_i 4 M_i/D_i). Equals 1 for
        two identical crystals; obeys reciprocity M1 t12 = M2 t21. Falls back to
        the k-integrated min(M1, M2) when k-resolved spectra are unavailable."""
        s1, s2 = self.s1, self.s2
        ok = (s1.M > 0) & (s2.M > 0) & (s1.D > 0) & (s2.D > 0)
        v1 = np.divide(4 * s1.M, s1.D, out=np.zeros_like(W), where=ok)
        v2 = np.divide(4 * s2.M, s2.D, out=np.zeros_like(W), where=ok)
        Z1, Z2 = s1.rho * v1, s2.rho * v2
        TZ = np.divide(4 * Z1 * Z2, (Z1 + Z2) ** 2, out=np.zeros_like(W), where=ok)
        Mmin = None
        if K_RESOLVED_SPECULAR and s1.crystal and s2.crystal:
            try:
                from phonons import channel_spectra
                M1k, V1k = channel_spectra(s1.name, s1.normal, nk=NK_SPEC, nz=NZ_SPEC, with_velocity=True)
                M2k, V2k = channel_spectra(s2.name, s2.normal, nk=NK_SPEC, nz=NZ_SPEC, with_velocity=True)
                if SMOOTH_THZ > 0:
                    M1k, M2k = _smooth(M1k, SMOOTH_THZ), _smooth(M2k, SMOOTH_THZ)
                # rescale so that each side integrates to its own full-mesh M (removes
                # the few-% quadrature difference between the two samplings)
                c1 = np.divide(s1.M, M1k.mean(axis=0), out=np.ones_like(W), where=M1k.mean(axis=0) > 0)
                c2 = np.divide(s2.M, M2k.mean(axis=0), out=np.ones_like(W), where=M2k.mean(axis=0) > 0)
                if K_IMPEDANCE:
                    # impedance factor per (k_par, w) from the group velocity of
                    # the crossing modes on each side (1D mode-matching estimate)
                    z1, z2 = s1.rho * V1k, s2.rho * V2k
                    tz = np.divide(4 * z1 * z2, (z1 + z2) ** 2, out=np.zeros_like(z1), where=(z1 + z2) > 0)
                    Mmin = (tz * np.minimum(M1k * c1, M2k * c2)).mean(axis=0)
                    TZ = np.ones_like(W)
                else:
                    Mmin = np.minimum(M1k * c1, M2k * c2).mean(axis=0)
                if LOWW_CONTINUUM:
                    # Long-wavelength limit: at low w only the few k_par columns near
                    # Gamma carry channels, and the finite k_par/q_z grid resolves
                    # them poorly (tau_spec -> 0 spuriously). In the continuum
                    # (acoustic) limit the transmission is frequency independent,
                    # so below the long-wavelength frequency w_b of the two sides
                    # (Side.w_debye) tau is held at its mean over [w_b, 1.5 w_b].
                    wb = max(s1.w_debye, s2.w_debye)
                    res = np.nonzero(W >= wb)[0]
                    if len(res) and wb > 0:
                        i0 = res[0]
                        win = (W >= wb) & (W <= 1.5 * wb) & ok
                        r0 = Mmin[win].sum() / s1.M[win].sum() if s1.M[win].sum() > 0 else 0.0
                        low = np.arange(len(W)) < i0
                        Mmin = np.where(low & ok, r0 * s1.M, Mmin)
            except Exception as exc:          # pragma: no cover
                import warnings
                warnings.warn(f"k-resolved specular transmission unavailable ({exc}); using k-integrated bound")
        if Mmin is None:
            Mmin = np.minimum(s1.M, s2.M)
        tau = np.divide(TZ * Mmin, s1.M, out=np.zeros_like(W), where=ok)
        if s1.name != s2.name:
            # channel counting cannot see eigenvector (polarisation/sublattice)
            # mismatch between two different crystals; the global efficiency is
            # calibrated against MACE-AGF on coherent pairs (calibrate_spec.py)
            tau = tau * spec_efficiency(s1.name, s2.name)
        return tau

    def specularity(self, w):
        """Fraction of phonons scattered specularly (parallel momentum
        conserved). A coherent interface with common in-plane period Lc
        conserves q_par modulo 2 pi/Lc; with a = smaller surface lattice
        parameter, a fraction p_inf = (a/Lc)^2 of the folded channels keeps the
        incident q_par (p_inf = 1 for a 1:1 pseudomorphic interface). Phonons
        with wavelength >> Lc see a smooth plane and are specular:
            p(w) = p_inf + (1 - p_inf) exp(-(Lc/lambda)^2),
        lambda = 2 pi v_T,min / w. Amorphous side -> Lc None -> p = 0."""
        if self.Lc is None:
            return np.zeros_like(w)
        p_inf = self.p_inf
        vs = min(self.s1.vT, self.s2.vT)
        lam = np.divide(2 * np.pi * vs, w, out=np.full_like(w, np.inf), where=w > 0)
        return p_inf + (1 - p_inf) * np.exp(-(self.Lc * 1e-10 / lam) ** 2)

    def bond_factor(self, w):
        if self.K is None:
            return np.ones_like(w)
        s1, s2 = self.s1, self.s2
        pn = np.interp(w, W, np.divide(s1.Mperp + s2.Mperp, s1.M + s2.M,
                                       out=np.full_like(W, 1 / 3), where=(s1.M + s2.M) > 0))
        if SPECTRAL_IMPEDANCE:
            # lattice impedance Z(w) = rho v_g(w), v_g = 4 M/D: vanishes at band
            # edges, where a soft interfacial spring no longer blocks transport
            def vg(s):
                return np.interp(w, W, np.divide(4 * s.M, s.D, out=np.zeros_like(W), where=s.D > 0))
            Z1, Z2 = s1.rho * vg(s1), s2.rho * vg(s2)
            ZL = ZT = (Z1, Z2)
        else:
            ZL = (s1.rho * s1.vL, s2.rho * s2.vL)
            ZT = (s1.rho * s1.vT, s2.rho * s2.vT)
        f = (pn * spring_factor(w, self.K[0], *ZL) / spring_factor(w, self.Kref[0], *ZL)
             + (1 - pn) * spring_factor(w, self.K[1], *ZT) / spring_factor(w, self.Kref[1], *ZT))
        return np.minimum(f, 1.0)

    def G_dmm_full(self, T):
        """Full-dispersion DMM (elastic, equilibrium temperatures, perfect bond)."""
        return np.sum(HBAR * W * dfdT(W, T) * self.s1.M * self.tau_dmm) * DW

    def G_elastic(self, T):
        return np.sum(HBAR * W * dfdT(W, T) * self.s1.M * self.tau_el * self.S) * DW

    def G_max(self, T):
        """Upper bound of the phonon conductance reachable with at most one
        three-phonon conversion step at the interface: the radiation limit of
        the softer side, or of the stiffer side restricted to w <= 2 w_c,soft
        (phonons that can split into two soft-side phonons), whichever is
        smaller; reduced by the flux-weighted interfacial bonding factor."""
        soft, stiff = (self.s1, self.s2) if self.s1.w_max <= self.s2.w_max else (self.s2, self.s1)
        wt = HBAR * W * dfdT(W, T)
        g_soft = np.sum(wt * soft.M) * DW
        g_stiff = np.sum(wt * stiff.M * (W <= 2 * soft.w_max)) * DW
        Sbar = np.sum(wt * soft.M * self.S) / np.sum(wt * soft.M)
        return min(g_soft, g_stiff) * Sbar

    def theta_soft(self):
        return HBAR * min(self.s1.w_max, self.s2.w_max) / KB

    def theta_stiff(self):
        return HBAR * max(self.s1.w_max, self.s2.w_max) / KB

    def P_inel(self, T, alpha, b):
        """Probability that an elastically reflected phonon is down/up-converted
        by interfacial anharmonic scattering into a transmissible mode before it
        leaves the interfacial region. Umklapp-type rate ~ T exp(-Theta/(b T))
        (Slack form), Theta = maximum phonon temperature of the STIFFER side:
        the conversion is the anharmonic decay of the stiff side's phonons
        that the soft side cannot accept elastically. (THETA_SIDE = "soft" is
        kept for comparison; it made soft metals such as Cu and Au strongly
        inelastic at 300 K and destroyed the Cu/TiC/diamond prediction.)"""
        if not self.sharp or alpha == 0:
            return 0.0
        if not (self.s1.crystal and self.s2.crystal) and not INELASTIC_AMORPHOUS:
            # Theta is the maximum phonon temperature of a computed crystal
            # spectrum; for an amorphous side it is only a fictitious Debye
            # cutoff, so the channel (observed/calibrated on sharp crystalline
            # interfaces) is not applied.
            return 0.0
        th = self.theta_stiff() if THETA_SIDE == "stiff" else self.theta_soft()
        return 1.0 - np.exp(-alpha * (T / th) * np.exp(-th / (b * T)))

    def G_pp(self, T, alpha=0.0, b=1.0, noneq=False):
        """Phonon-phonon conductance. Elastic channel (specular/diffuse mix of
        full-dispersion transmissions, times the bonding factor) plus the
        inelastic conversion channel:
            G = G_el + P(T) (G_max - G_el)
        noneq=True applies the Zeng-Chen local-equilibrium correction
        (not used by the frozen model; kept for comparison)."""
        Gel = self.G_elastic(T)
        GL = Gel + self.P_inel(T, alpha, b) * max(self.G_max(T) - Gel, 0.0)
        if not noneq:
            return GL
        R = 1 / GL - 0.5 * (1 / self.s1.G_rad(T) + 1 / self.s2.G_rad(T))
        return 1 / R if R > 1e-3 / GL else np.inf          # identical sides: no resistance




def kappa_min(side, T, n_int=400):
    """Cahill-Watson-Pohl minimum thermal conductivity of an amorphous solid
    (PRB 46, 6131 (1992)) from its sound velocities and atom density:
        k = (pi/6)^(1/3) kB n^(2/3) sum_i v_i (T/Th_i)^2 int_0^{Th_i/T} x^3 e^x/(e^x-1)^2 dx,
        Th_i = v_i (hbar/kB) (6 pi^2 n)^(1/3).
    Used for the conduction of thin glassy films (native/thermal oxide) at any T."""
    n = side.n_atoms
    k = 0.0
    for v in (side.vL, side.vT, side.vT):
        th = v * HBAR / KB * (6 * np.pi ** 2 * n) ** (1 / 3)
        xm = th / T
        x = np.linspace(1e-6, xm, n_int)
        f = x ** 3 * np.exp(x) / np.expm1(x) ** 2
        k += v * (T / th) ** 2 * np.trapezoid(f, x)
    return (np.pi / 6) ** (1 / 3) * KB * n ** (2 / 3) * k
