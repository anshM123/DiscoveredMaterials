"""Thermal boundary conductance: baselines (Debye AMM / DMM) and the
transferable model described in README.md.

All conductances are returned in W m^-2 K^-1 unless stated otherwise.
Angular frequency grid W (rad/s) is shared with phonons.py.
"""
import numpy as np
from phonons import W, DW, load, T_GRID
from materials import CRYSTALS, AMORPHOUS, METALS

HBAR = 1.054571817e-34
KB = 1.380649e-23
AMU = 1.66053906660e-27
THETA_SIDE = "stiff"
INELASTIC_AMORPHOUS = False   # inelastic channel only for crystal/crystal interfaces


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
            # in-plane surface lattice parameter for (111)/(0001) normals (Angstrom)
            self.a_surf = float(np.linalg.norm(d["a_cell"][0]))
            self.rho = float(d["rho"])
            self.vL, self.vT = float(d["vL"]), float(d["vT"])
            self.n_atoms = float(d["n_atoms"])
            self.m_avg = float(d["m_avg"])
            self.crystal = True
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
            # isotropic solid: flux-weighted fraction of displacement along n is 1/2
            self.Mperp = 0.5 * ML + 0.5 * MT
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
    tau = np.where(ok, 4 * Z1 * Z2 * c1 * c2 / (Z1 * c1 + Z2 * c2) ** 2, 0.0)
    return 2 * np.sum(tau * c1 * np.sin(th)) * (np.pi / 2) / n


def G_AMM_debye(s1, s2, T):
    """Acoustic mismatch model, Debye dispersion (Little 1959), evaluated from side 1."""
    G = 0.0
    for (v1, wc1, g), (v2, wc2, _) in zip(s1.debye_branches(), s2.debye_branches()):
        a = _amm_angle_avg(s1.rho * v1, s2.rho * v2, v1, v2)
        M1 = np.where(W <= wc1, g * W ** 2 / (2 * np.pi ** 2 * v1 ** 3) * v1 / 4, 0.0)
        G += np.sum(HBAR * W * dfdT(W, T) * M1) * DW * a
    return G


# ---------------------------------------------------------------- the model
COHERENCE_STRAIN = 0.03   # max misfit strain for a pseudomorphic (coherent) match


def coincidence_length(s1, s2, nmax=12, tol=COHERENCE_STRAIN):
    """Smallest common in-plane period (Angstrom) of two hexagonal surface
    lattices in the aligned orientation (close-packed rows parallel, the
    usual epitaxial relation for (111)/(0001) stacks) within misfit strain
    tol; None if a side is amorphous."""
    if s1.a_surf is None or s2.a_surf is None:
        return None
    best = None
    for f1 in (1.0,):
        for f2 in (1.0,):
            a1, a2 = s1.a_surf * f1, s2.a_surf * f2
            for n1 in range(1, nmax + 1):
                for n2 in range(1, nmax + 1):
                    if abs(n1 * a1 / (n2 * a2) - 1) <= tol:
                        L = max(n1 * a1, n2 * a2)
                        if best is None or L < best:
                            best = L
    return best


def spring_factor(w, K, Z1, Z2):
    """Energy transmission of a plane wave through an interfacial spring K
    (N m^-3) relative to a welded contact."""
    zr = Z1 * Z2 / (Z1 + Z2)
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

    def __init__(self, s1, s2, K=None, Kref=None, sharp=True, Lc="auto"):
        self.s1, self.s2 = s1, s2
        self.K, self.Kref, self.sharp = K, Kref, sharp
        M1, M2 = s1.M, s2.M
        self.tau_dmm = np.divide(M2, M1 + M2, out=np.zeros_like(W), where=(M1 + M2) > 0)
        self.Lc = coincidence_length(s1, s2) if Lc == "auto" else Lc
        self.tau_spec = self.specular_transmission()
        self.p_spec = self.specularity(W)
        self.tau_el = self.p_spec * self.tau_spec + (1 - self.p_spec) * self.tau_dmm
        self.S = self.bond_factor(W)

    def specular_transmission(self):
        """Momentum-conserving (specular) limit: impedance factor times the
        fraction of channels available on the other side. Equals 1 for two
        identical crystals; obeys reciprocity M1 t12 = M2 t21."""
        s1, s2 = self.s1, self.s2
        ok = (s1.M > 0) & (s2.M > 0) & (s1.D > 0) & (s2.D > 0)
        v1 = np.divide(4 * s1.M, s1.D, out=np.zeros_like(W), where=ok)
        v2 = np.divide(4 * s2.M, s2.D, out=np.zeros_like(W), where=ok)
        Z1, Z2 = s1.rho * v1, s2.rho * v2
        TZ = np.divide(4 * Z1 * Z2, (Z1 + Z2) ** 2, out=np.zeros_like(W), where=ok)
        return np.divide(TZ * np.minimum(s1.M, s2.M), s1.M, out=np.zeros_like(W), where=ok)

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
        a = min(self.s1.a_surf, self.s2.a_surf)
        p_inf = min(1.0, (a / self.Lc) ** 2)
        vs = min(self.s1.vT, self.s2.vT)
        lam = np.divide(2 * np.pi * vs, w, out=np.full_like(w, np.inf), where=w > 0)
        return p_inf + (1 - p_inf) * np.exp(-(self.Lc * 1e-10 / lam) ** 2)

    def bond_factor(self, w):
        if self.K is None:
            return np.ones_like(w)
        s1, s2 = self.s1, self.s2
        pn = np.interp(w, W, np.divide(s1.Mperp + s2.Mperp, s1.M + s2.M,
                                       out=np.full_like(W, 1 / 3), where=(s1.M + s2.M) > 0))
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
        return 1 / R if R > 0 else np.inf


