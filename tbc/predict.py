"""High-level predictions for the three reference systems and arbitrary
stacks. Global parameters live in PARAMS and are frozen after calibration
(see calibrate.py); nothing is fitted per specimen.
"""
import json, os
import numpy as np
from model import Side, Interface, G_DMM_debye, G_AMM_debye
from materials import METALS, FILM_K

HERE = os.path.dirname(os.path.abspath(__file__))
PARAM_FILE = os.path.join(HERE, "params.json")


def load_params():
    with open(PARAM_FILE) as f:
        return json.load(f)


def G_ep(metal, T, d=None):
    """Electron-phonon conductance (W m^-2 K^-1).
    Semi-infinite metal: Majumdar-Reddy sqrt(g kp(T)).
    Thin layer of thickness d (m) << e-ph length: volumetric coupling g*d."""
    p = metal_params(metal)
    if p is None:
        import warnings
        warnings.warn(f"No electron-phonon data for metal '{metal}' in materials.METALS: "
                      "the electron-phonon resistance is omitted (G is then an upper bound). "
                      "Add g and kp for this metal to include it.", stacklevel=2)
        return np.inf
    if d is not None:
        return p["g"] * d
    return np.sqrt(p["g"] * kp_metal(metal, T))


def g_allen(metal):
    """Allen's estimate of g (W m^-3 K^-1) from the Sommerfeld coefficient, the
    electron-phonon coupling constant and the MACE phonon DOS (materials.py)."""
    from materials import SOMMERFELD, LAMBDA_EP
    from phonons import load, W
    hbar, kB, NA = 1.054571817e-34, 1.380649e-23, 6.02214076e23
    d = load(metal)
    s = W > 0
    w2 = np.sum(d["D"][s] * W[s]) / np.sum(d["D"][s] / W[s])
    lam = LAMBDA_EP[metal]
    gamma_b = SOMMERFELD[metal] * 1e-3 / (d["V"] / d["nat"] * NA) / (1 + lam)
    return 3 * hbar * gamma_b * lam * w2 / (np.pi * kB)


_METAL_CACHE = {}


def metal_params(metal):
    """g, kp, kp_pp for a metal: measured g (materials.METALS) where available,
    else Allen's estimate (with a warning); None if nothing is known."""
    from materials import SOMMERFELD, LAMBDA_EP, KP_EXTRA, KP_NOMINAL, CRYSTALS
    if metal in METALS:
        return METALS[metal]
    if metal in _METAL_CACHE:
        return _METAL_CACHE[metal]
    if metal in SOMMERFELD and metal in LAMBDA_EP and metal in CRYSTALS:
        import warnings
        kp, kp_pp = KP_EXTRA.get(metal, KP_NOMINAL)
        p = dict(g=g_allen(metal), kp=kp, kp_pp=kp_pp, estimated=True)
        warnings.warn(f"{metal}: no measured electron-phonon coupling; Allen estimate "
                      f"g = {p['g']:.2e} W/m^3K (factor-2 uncertainty)", stacklevel=3)
        _METAL_CACHE[metal] = p
        return p
    return None


def kp_metal(metal, T):
    """Lattice thermal conductivity of a metal vs T. Phonon-phonon part ~1/T,
    phonon-electron part T-independent (Matthiessen), both at 300 K from
    Tong et al. PRB 100, 144306 (2019), Table 2 (arXiv:1907.01364)."""
    p = metal_params(metal)
    k_pe = 1.0 / (1.0 / p["kp"] - 1.0 / p["kp_pp"])
    return 1.0 / (T / 300.0 / p["kp_pp"] + 1.0 / k_pe)


def interface(a, b, na=None, nb=None, sharp=True, bonding=None, epitaxial=True):
    s1, s2 = Side(a, na), Side(b, nb)
    K = Kref = None
    if bonding is not None:
        K, Kref = bonding
    return Interface(s1, s2, K=K, Kref=Kref, sharp=sharp, epitaxial=epitaxial)


def G_metal_nonmetal(metal, other, T, P, nm=None, no=None, sharp=True, bonding=None):
    """Metal on a non-metal: e-ph resistance in series with phonon interface."""
    I = interface(metal, other, nm, no, sharp=sharp, bonding=bonding)
    Gpp = I.G_pp(T, alpha=P["alpha"], b=P["b"], noneq=P.get("noneq", False))
    return 1.0 / (1.0 / G_ep(metal, T) + 1.0 / Gpp), Gpp


def baselines(a, b, T, na=None, nb=None):
    s1, s2 = Side(a, na), Side(b, nb)
    return dict(AMM=G_AMM_debye(s1, s2, T), DMM=G_DMM_debye(s1, s2, T))


# ------------------------------------------------------------------ stacks
def sandell_stack(adhesion, T=300.0, P=None, bond=None, t_ox=2.05e-9, t_ad=None):
    """Combined boundary conductance of Au/[adhesion]/PMMA/SiO2(native)/Si, as
    extracted by Sandell et al.: everything except the PMMA bulk resistance.
    Electrons in Au deliver heat to the lattice of Au (Majumdar-Reddy) and, in
    parallel, to the lattice of the thin adhesion layer (volumetric g*d)."""
    P = P or load_params()
    bond = bond or {}
    Gep = G_ep("Au", T)
    top = "Au"
    R = 0.0
    if adhesion is not None:
        t_ad = t_ad or {"Ti": 2.4e-9, "Ni": 2.0e-9}[adhesion]
        Gep += G_ep(adhesion, T, d=t_ad)
        top = adhesion
        R += t_ad / FILM_K[adhesion]
    R += 1 / Gep
    parts = {}
    for a, b in [(top, "PMMA"), ("PMMA", "SiO2"), ("SiO2", "Si")]:
        I = interface(a, b, bonding=bond.get((a, b)))
        parts[f"{a}/{b}"] = I.G_pp(T, P["alpha"], P["b"])
        R += 1 / parts[f"{a}/{b}"]
    R += t_ox / film_k("SiO2", T)
    parts["G_ep"] = Gep
    return 1 / R, parts


def chang_ratio(T=300.0, P=None, bond=None, t_tic=10e-9):
    """G(Cu/TiC(10nm)/diamond) / G(Cu/diamond). Heat path: Cu electrons ->
    Cu lattice (Majumdar-Reddy) -> phonon interfaces and TiC film in series."""
    P = P or load_params()
    bond = bond or {}
    Gep = G_ep("Cu", T)
    # magnetron-sputtered, polycrystalline Cu and TiC films (Chang et al.):
    # non-epitaxial -> diffuse scattering
    g = lambda a, b: interface(a, b, "111", "111", bonding=bond.get((a, b)), epitaxial=False).G_pp(T, P["alpha"], P["b"])
    G0 = 1 / (1 / Gep + 1 / g("Cu", "C"))
    G1 = 1 / (1 / Gep + 1 / g("Cu", "TiC") + t_tic / FILM_K["TiC"] + 1 / g("TiC", "C"))
    return G1 / G0, dict(G_Cu_C=G0, G_Cu_TiC_C=G1, Cu_C=g("Cu", "C"), Cu_TiC=g("Cu", "TiC"), TiC_C=g("TiC", "C"), G_ep=Gep)


# ------------------------------------------------------------------ bonding
BOND_KEYS = {("Al", "Si"): "Al_Si", ("Al", "GaN"): "Al_GaN", ("Cu", "TiC"): "Cu_TiC",
             ("TiC", "C"): "TiC_C", ("Au", "PMMA"): "Au_PMMA", ("Ti", "PMMA"): "Ti_PMMA",
             ("Ni", "PMMA"): "Ni_PMMA", ("PMMA", "SiO2"): "PMMA_SiO2",
             ("Ag", "PMMA"): "Ag_PMMA", ("Pt", "PMMA"): "Pt_PMMA"}


MOLECULAR_CONTACTS_IDEAL = True


def bonding_for(a, b, na=None, nb=None, key=None):
    """(K, Kref) for an interface from the MACE bonding cache.
    The bond-strength ratio r = K_n / K_ref,n (normal rigid-shift stiffness)
    is applied to both polarisations: K = r * K_ref. The uniform-shear
    stiffness is not used because it averages to ~0 for incommensurate
    interfaces and does not represent dynamic coupling of short phonons.
    Returns None (ideal bond) if the interface has not been computed."""
    from bonding import load_cache, Kref
    from materials import AMORPHOUS, physisorbs
    c = load_cache()
    key = key or BOND_KEYS.get((a, b)) or BOND_KEYS.get((b, a))
    if MOLECULAR_CONTACTS_IDEAL and ("PMMA" in (a, b)):
        # molecule/surface contacts. Chemisorbed contacts (reactive metals, the
        # hydroxylated oxide surface): MACE-OMAT-0 (inorganic training data) does not
        # describe them (it fails to bind the PMMA carbonyl to Ni, which XPS shows),
        # so they are treated as ideally bonded. Physisorbed contacts (noble metals,
        # materials.physisorbs): dispersion-only binding, which MACE + D3 does
        # describe -> its rigid-shift stiffness is used.
        other = b if a == "PMMA" else a
        if other in AMORPHOUS or not physisorbs(other):
            return None
        kr = Kref(a, b, na, nb)
        if key is None or key not in c:
            # no MACE scan for this metal: dispersion-bonded contact with the
            # bond-strength ratio of Au/PMMA (closest computed analogue)
            import warnings
            warnings.warn(f"{a}/{b}: physisorbed contact without a MACE scan; using the Au/PMMA "
                          "bond-strength ratio (run bonding.py for this pair)", stacklevel=2)
            r = c["Au_PMMA"]["K_n"] / Kref("Au", "PMMA")["K_n"]
            return ((r * kr["K_n"], r * kr["K_t"]), (kr["K_n"], kr["K_t"]))
    if key is None or key not in c:
        return None
    kr = Kref(a, b, na, nb)
    r = c[key]["K_n"] / kr["K_n"]
    return ((r * kr["K_n"], r * kr["K_t"]), (kr["K_n"], kr["K_t"]))


def film_k(name, T):
    """Thermal conductivity of a thin film. Glassy films: measured 300 K value
    (materials.FILM_K) times the temperature dependence of the Cahill-Watson-
    Pohl minimum conductivity, k(T) = k(300) k_min(T)/k_min(300); crystalline
    films: constant."""
    from model import kappa_min
    from materials import AMORPHOUS
    k300 = FILM_K[name]
    if name in AMORPHOUS:
        s = Side(name)
        return k300 * kappa_min(s, T) / kappa_min(s, 300.0)
    return k300
