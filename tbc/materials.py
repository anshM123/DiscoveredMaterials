"""Material definitions.

Crystalline materials are described only by an initial structure guess; the
lattice is relaxed and all phonon properties are computed with the MACE-OMAT-0
universal interatomic potential (see phonons.py). Nothing here is fitted to
thermal-boundary-conductance data.

The few non-phonon inputs (electron-phonon coupling, phonon thermal
conductivity of metals, Born charges, properties of amorphous polymers) are
literature values; each carries its source.
"""
from dataclasses import dataclass, field
import re
import numpy as np
from ase.build import bulk
from ase import Atoms


def wurtzite(sym_a, sym_b, a, c, u):
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
    sp = [[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, 0.5],
          [1 / 3, 2 / 3, u], [2 / 3, 1 / 3, 0.5 + u]]
    return Atoms([sym_a, sym_a, sym_b, sym_b], cell=cell, scaled_positions=sp, pbc=True)


# name -> (structure builder, phonopy supercell, q-mesh, interface normals
#          available as {label: cartesian direction in the conventional frame})
CUBIC_111 = {"111": [1, 1, 1], "100": [0, 0, 1], "110": [1, 1, 0]}
HEX_0001 = {"0001": [0, 0, 1]}

def corundum(sym_m="Al", a=4.759, c=12.991, zm=0.3522, xo=0.3064):
    """alpha-Al2O3 (R-3c), hexagonal setting (30 atoms), c along z."""
    from ase.spacegroup import crystal
    return crystal([sym_m, "O"], basis=[(0, 0, zm), (xo, 0, 0.25)], spacegroup=167,
                   cellpar=[a, a, c, 90, 90, 120], primitive_cell=False)


def sic_4h(a=3.081, c=10.085):
    """4H-SiC (P6_3mc), stacking ABCB, c along z; C above Si by 3c/16."""
    from ase import Atoms
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
    st = [(0, 0, 0.0), (1 / 3, 2 / 3, 0.25), (2 / 3, 1 / 3, 0.5), (1 / 3, 2 / 3, 0.75)]
    sp = [(x, y, z) for x, y, z in st] + [(x, y, z + 3 / 16) for x, y, z in st]
    return Atoms(["Si"] * 4 + ["C"] * 4, cell=cell, scaled_positions=sp, pbc=True)


def beta_si3n4(a=7.608, c=2.911):
    """beta-Si3N4 (P6_3/m), c along z."""
    from ase.spacegroup import crystal
    return crystal(["Si", "N", "N"], basis=[(0.1733, 0.7694, 0.25), (1 / 3, 2 / 3, 0.25), (0.3323, 0.0301, 0.25)],
                   spacegroup=176, cellpar=[a, a, c, 90, 90, 120], primitive_cell=False)


CRYSTALS = {
    "Al":  dict(build=lambda: bulk("Al", "fcc", a=4.04), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Cu":  dict(build=lambda: bulk("Cu", "fcc", a=3.61), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Au":  dict(build=lambda: bulk("Au", "fcc", a=4.08), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Ni":  dict(build=lambda: bulk("Ni", "fcc", a=3.52), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Ti":  dict(build=lambda: bulk("Ti", "hcp", a=2.95, c=4.68), sc=[4, 4, 3], mesh=[32, 32, 20], normals=HEX_0001),
    "Si":  dict(build=lambda: bulk("Si", "diamond", a=5.43), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111),
    "Ge":  dict(build=lambda: bulk("Ge", "diamond", a=5.658), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111),
    "C":   dict(build=lambda: bulk("C", "diamond", a=3.567), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111),
    "TiC": dict(build=lambda: bulk("TiC", "rocksalt", a=4.33), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111),
    "GaN": dict(build=lambda: wurtzite("Ga", "N", 3.19, 5.19, 0.377), sc=[4, 4, 3], mesh=[28, 28, 18], normals=HEX_0001,
                # Born effective charge / eps_inf for LO-TO splitting (isotropic approx.)
                # Z*=2.65, eps_inf=5.35: Bernardini & Fiorentini PRB 56, R10024 (1997); Davydov PRB 58, 12899 (1998)
                nac=dict(Z=2.65, eps=5.35)),
    # Born charges / eps_inf (isotropic approximation):
    #   SiC Z*=2.7, eps=6.5 (Karch et al., PRB 50, 17054 (1994)); GaAs Z*=2.07, eps=10.9
    #   and AlAs Z*=2.17, eps=8.2 (Giannozzi et al., PRB 43, 7231 (1991))
    "SiC4H": dict(build=lambda: sic_4h(), sc=[4, 4, 2], mesh=[24, 24, 8], normals=HEX_0001,
                  nac=dict(Z=2.7, eps=6.5, cation="Si")),
    "SiC3C": dict(build=lambda: bulk("SiC", "zincblende", a=4.36), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                  nac=dict(Z=2.7, eps=6.5, cation="Si")),
    "GaAs": dict(build=lambda: bulk("GaAs", "zincblende", a=5.653), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                 nac=dict(Z=2.07, eps=10.9, cation="Ga")),
    "AlAs": dict(build=lambda: bulk("AlAs", "zincblende", a=5.661), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                 nac=dict(Z=2.17, eps=8.2, cation="Al")),
    # ZnO wurtzite: Z*=2.1, eps_inf=3.7 (Zaoui et al. / standard values; isotropic approx.)
    "ZnO": dict(build=lambda: wurtzite("Zn", "O", 3.25, 5.207, 0.382), sc=[4, 4, 3], mesh=[28, 28, 18], normals=HEX_0001,
                nac=dict(Z=2.1, eps=3.7, cation="Zn")),
    # ---- additional microelectronics materials (structure guesses; MACE-relaxed) ----
    "TiN": dict(build=lambda: bulk("TiN", "rocksalt", a=4.24), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111),
    "MgO": dict(build=lambda: bulk("MgO", "rocksalt", a=4.21), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                nac=dict(Z=1.98, eps=2.95, cation="Mg")),
    "W":   dict(build=lambda: bulk("W", "bcc", a=3.165), sc=[5, 5, 5], mesh=[40, 40, 40], normals=CUBIC_111),
    "Mo":  dict(build=lambda: bulk("Mo", "bcc", a=3.147), sc=[5, 5, 5], mesh=[40, 40, 40], normals=CUBIC_111),
    "Ta":  dict(build=lambda: bulk("Ta", "bcc", a=3.301), sc=[5, 5, 5], mesh=[40, 40, 40], normals=CUBIC_111),
    "Pt":  dict(build=lambda: bulk("Pt", "fcc", a=3.92), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Ag":  dict(build=lambda: bulk("Ag", "fcc", a=4.086), sc=[4, 4, 4], mesh=[40, 40, 40], normals=CUBIC_111),
    "Co":  dict(build=lambda: bulk("Co", "hcp", a=2.507, c=4.07), sc=[4, 4, 3], mesh=[32, 32, 20], normals=HEX_0001),
    "Ru":  dict(build=lambda: bulk("Ru", "hcp", a=2.706, c=4.282), sc=[4, 4, 3], mesh=[32, 32, 20], normals=HEX_0001),
    "InP": dict(build=lambda: bulk("InP", "zincblende", a=5.869), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                nac=dict(Z=2.55, eps=9.6, cation="In")),
    "GaP": dict(build=lambda: bulk("GaP", "zincblende", a=5.451), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                nac=dict(Z=2.15, eps=9.1, cation="Ga")),
    "cBN": dict(build=lambda: bulk("BN", "zincblende", a=3.615), sc=[4, 4, 4], mesh=[36, 36, 36], normals=CUBIC_111,
                nac=dict(Z=1.9, eps=4.5, cation="B")),
    "Si3N4": dict(build=lambda: beta_si3n4(), sc=[2, 2, 4], mesh=[12, 12, 24], normals=HEX_0001),
    # sapphire; Born charge Z(Al)=2.9, eps_inf=3.1: Heid, Strauch & Bohnen, PRB 61, 12059 (2000)
    "Al2O3": dict(build=lambda: corundum(), sc=[2, 2, 1], mesh=[14, 14, 5], normals=HEX_0001,
                  nac=dict(Z=2.9, eps=3.1, cation="Al")),
    "AlN": dict(build=lambda: wurtzite("Al", "N", 3.11, 4.98, 0.382), sc=[4, 4, 3], mesh=[28, 28, 18], normals=HEX_0001,
                nac=dict(Z=2.65, eps=4.4)),
}

# Electron-phonon coupling factor g (W m^-3 K^-1) and lattice thermal
# conductivity of metals, used for the Majumdar-Reddy electron-phonon
# conductance G_ep = sqrt(g * kp(T)).
#   kp   : phonon thermal conductivity at 300 K with phonon-phonon AND
#          phonon-electron scattering; kp_pp: phonon-phonon only. Both DFT,
#          Tong, Li, Ruan & Bao, PRB 100, 144306 (2019), Table 2 (arXiv:1907.01364).
#   g    : room-temperature experimental values (W m^-3 K^-1):
#          Al 2.45e17, Au 2.9e16 : Hostetler, Smith & Norris, Appl. Opt. 38, 3614 (1999), Table 1
#          Cu 1.0e17             : Elsayed-Ali et al., PRL 58, 1212 (1987)
#          Ni 8.6e17             : Kang & Choi, J. Magn. Magn. Mater. 514, 167156 (2020)
#          Ti 1.3e18             : room-T thermoreflectance value quoted by
#                                  Lin, Zhigilei & Celli, PRB 77, 075133 (2008)
#          Note: for Al the equilibrium (DFT) g is ~2x the pump-probe value
#          (Brown et al., PRB 94, 075120 (2016)); g carries +-50% uncertainty.
METALS = {
    "Al": dict(g=2.45e17, kp=8.14, kp_pp=10.67),
    "Cu": dict(g=1.0e17, kp=17.61, kp_pp=19.64),
    "Au": dict(g=2.9e16, kp=2.72, kp_pp=3.05),
    "Ni": dict(g=8.6e17, kp=14.50, kp_pp=28.92),
    "Ti": dict(g=1.3e18, kp=6.68, kp_pp=11.92),
}

# Metals without a measured g: Allen's formula (PRL 59, 1460 (1987)),
#   g = 3 hbar gamma_b lambda <w^2> / (pi kB),  gamma_b = gamma / (1 + lambda),
# with <w^2> = int F w / int F/w from the MACE phonon DOS F (constant alpha^2).
# Checked on the five metals with measured g above: ratio 0.42 (Au) .. 2.7 (Ni),
# geometric mean 0.91 -> a factor-2 estimate (predict.g_allen).
#   gamma  : electronic heat-capacity coefficient, mJ mol^-1 K^-2 (Kittel,
#            Introduction to Solid State Physics, 8th ed., Table 6.2)
#   lambda : electron-phonon coupling constant. Tong et al. (2019) Table 2
#            (experimental column; Co: DFT); W, Mo, Ta, Ru: McMillan, Phys. Rev.
#            167, 331 (1968), Table IV
SOMMERFELD = {"Al": 1.35, "Cu": 0.695, "Ag": 0.646, "Au": 0.729, "Ni": 7.02, "Ti": 3.35,
              "Pt": 6.8, "Co": 4.73, "W": 1.3, "Mo": 2.0, "Ta": 5.9, "Ru": 3.3}
LAMBDA_EP = {"Al": 0.43, "Cu": 0.14, "Ag": 0.13, "Au": 0.15, "Ni": 0.31, "Ti": 0.38,
             "Pt": 0.66, "Co": 0.42, "W": 0.28, "Mo": 0.41, "Ta": 0.65, "Ru": 0.38}
# lattice conductivity (kp, kp_pp) of further metals, Tong et al. (2019) Table 2;
# metals absent there get a nominal 10 / 15 W/mK (Tong: 2-18 W/mK for elemental
# metals). G_ep ~ sqrt(kp), and for these strongly coupled metals G_ep >> G_phonon.
KP_EXTRA = {"Ag": (5.59, 6.03), "Pt": (6.28, 8.72), "Co": (12.20, 20.60)}
KP_NOMINAL = (10.0, 15.0)


# Free-electron Fermi velocities (m/s), Ashcroft & Mermin, Solid State Physics,
# Table 2.1, for the electron diffuse-mismatch channel of metal/metal
# interfaces (predict.G_edmm). d-band metals without a tabulated free-electron
# value get a nominal 1.0e6 m/s (flagged).
FERMI_V = {"Al": 2.03e6, "Cu": 1.57e6, "Ag": 1.39e6, "Au": 1.40e6}
FERMI_V_NOMINAL = 1.0e6


# Standard formation enthalpy of each metal's most stable oxide per mole of O
# (kJ/mol O; CRC Handbook / NIST-JANAF). Used to classify metal / molecular-
# solid contacts: metals whose oxide is barely stable (> -100 kJ/mol O: Au, Ag,
# Pt, Pd) do not form metal-O bonds with ester/carbonyl oxygen and contact a
# polymer by dispersion only (physisorption); all others chemisorb (XPS
# reactivity series for metals on PMMA: Ti > Cr > Al > Ni > Cu >> Ag, Au).
OXIDE_DHF_PER_O = {"Au": -1.0, "Ag": -31.1, "Pt": -67.0, "Pd": -85.4, "Cu": -168.6,
                   "Ni": -239.7, "Co": -237.9, "Ru": -152.5, "W": -281.0, "Mo": -248.4,
                   "Ta": -409.2, "Ti": -472.0, "Al": -558.6, "Cr": -379.9, "Fe": -274.7}
PHYSISORB_DHF = -100.0


def physisorbs(metal):
    """True if a molecular solid binds to this metal by dispersion only."""
    return OXIDE_DHF_PER_O.get(metal, -1e9) > PHYSISORB_DHF


# Thin-film thermal conductivities used for series (bulk) resistances of
# interlayers. TiC: ~20 W/mK (bulk TiC 21-33 W/mK); PMMA: 0.20 W/mK (accepted
# bulk value, Sandell et al. 2020); SiO2: 1.38 W/mK (fused silica).
FILM_K = {"TiC": 20.0, "PMMA": 0.20, "SiO2": 1.38,
          # 2 nm adhesion films: values used by Sandell et al. (2020), Table 1
          "Ti": 8.2, "Ni": 52.0}
# Fallback for other crystalline films: room-temperature BULK conductivity
# (W/mK, CRC Handbook), an upper bound for nm films, whose bulk resistance is
# in any case small next to their interfaces. Amorphous films without a
# measured value use the Cahill-Watson-Pohl minimum conductivity.
BULK_K = {"Al": 237, "Cu": 401, "Ag": 429, "Au": 317, "Pt": 72, "Ni": 91, "Co": 100, "Ti": 22, "W": 173,
          "Mo": 138, "Ta": 57, "Ru": 117, "TiN": 29, "TiC": 21, "Si": 148, "Ge": 60, "C": 2000, "GaN": 230,
          "AlN": 285, "Al2O3": 35, "MgO": 55, "SiC4H": 370, "SiC3C": 360, "GaAs": 55, "AlAs": 91, "InP": 68,
          "GaP": 110, "ZnO": 50, "cBN": 740, "Si3N4": 30}


@dataclass
class Amorphous:
    """Debye-like description of an amorphous solid / glassy polymer."""
    name: str
    rho: float          # kg/m^3
    vL: float           # m/s
    vT: float           # m/s
    n_atoms: float      # atoms per m^3
    m_avg: float        # amu, mean atomic mass
    omega_max_THz: float  # highest frequency that carries heat (THz, cycles)
    notes: str = ""


# PMMA: rho = 1180 kg/m3, vL = 2690 m/s, vT = 1340 m/s (ultrasonic, room T).
# Atom density: C5O2H8 repeat unit (100.1 g/mol, 15 atoms).
# Debye cutoffs of amorphous solids are taken from the atomic number density
# (Swartz & Pohl convention), the same rule for every amorphous material.
_pmma_units = 1180 / (100.12e-3 / 6.02214076e23)
PMMA = Amorphous("PMMA", rho=1180.0, vL=2690.0, vT=1340.0,
                 n_atoms=_pmma_units * 15, m_avg=100.12 / 15,
                 omega_max_THz=0.0, notes="Debye cutoff from atom density")

# Amorphous SiO2 (native oxide / thermal oxide): rho = 2200 kg/m3,
# vL = 5953 m/s, vT = 3743 m/s (fused silica, ultrasonic); SiO2 = 3 atoms / 60.08 g/mol.
_sio2_units = 2200 / (60.08e-3 / 6.02214076e23)
SIO2 = Amorphous("SiO2", rho=2200.0, vL=5953.0, vT=3743.0,
                 n_atoms=_sio2_units * 3, m_avg=60.08 / 3,
                 omega_max_THz=0.0, notes="Debye cutoff from atom density")
# Amorphous Al2O3 (ALD / anodic alumina, the common gate and passivation
# dielectric): rho = 3050 kg/m3; sound velocities from the elastic constants
# typical of ALD alumina, E ~ 170 GPa and Poisson ratio ~0.24 (nanoindentation
# range E = 150-180 GPa): vL = 8100, vT = 4740 m/s. (Cross-check: the MACE
# melt-quenched glass of glass.py gives 7810 / 5240 m/s, i.e. -4% / +10%; for
# a-SiO2 it gives 6910 / 4750 vs 5953 / 3743 m/s measured for fused silica.)
_al2o3_units = 3050 / (101.961e-3 / 6.02214076e23)
AAL2O3 = Amorphous("aAl2O3", rho=3050.0, vL=8100.0, vT=4740.0,
                   n_atoms=_al2o3_units * 5, m_avg=101.961 / 5,
                   omega_max_THz=0.0, notes="Debye cutoff from atom density")
AMORPHOUS = {"PMMA": PMMA, "SiO2": SIO2, "aAl2O3": AAL2O3}


def register_structure(path, miller, name=None):
    """Add a crystal from a structure file (CIF, POSCAR, extxyz, ... via ASE)
    to CRYSTALS, with interface normal given by Miller indices (hkl) of the
    file's cell (or 'hkil' / '0001' style for hexagonal cells). Supercell and
    q-mesh are chosen from the cell size. No Born charges are set: LO-TO
    splitting of polar crystals is omitted (warned)."""
    import os, hashlib
    from ase.io import read
    from ase import Atoms
    atoms0 = read(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    if name is None:
        # never shadow a library material or a different structure with the same
        # file name: the cache key carries a hash of the structure
        h = hashlib.sha1(np.round(np.concatenate([atoms0.cell.array.ravel(), atoms0.get_scaled_positions().ravel()]), 5)
                         .tobytes() + "".join(atoms0.get_chemical_symbols()).encode()).hexdigest()[:8]
        name = f"{stem}_{h}"
    # reduce to the primitive cell (keeps the Cartesian orientation), so that the
    # interface-plane lattice and k_par matching use the true plane Bravais lattice
    try:
        import spglib
        cell = (atoms0.cell.array, atoms0.get_scaled_positions(), atoms0.get_atomic_numbers())
        lat, pos, num = spglib.standardize_cell(cell, to_primitive=True, no_idealize=True, symprec=1e-3)
        atoms = Atoms(numbers=num, cell=lat, scaled_positions=pos, pbc=True)
    except Exception:
        atoms = atoms0
    m = str(miller)
    # Miller indices: digits with optional minus signs, e.g. 111, 1-10, -1-12, 0001
    toks = re.findall(r"-?\d", m) if re.fullmatch(r"(-?\d)+", m) else None
    digits = [int(t) for t in toks] if toks else None
    if digits is None or len(digits) not in (3, 4):
        raise ValueError(f"Miller indices '{miller}' not understood (use e.g. 001, 111, 0001)")
    if len(digits) == 4:
        digits = [digits[0], digits[1], digits[3]]
    rec = np.linalg.inv(atoms0.cell.array).T          # Miller indices refer to the file's cell
    nvec = np.asarray(digits, float) @ rec
    lengths = atoms.cell.lengths()
    sc = [int(np.ceil(10.0 / L)) for L in lengths]
    mesh = [int(np.clip(np.ceil(120.0 / L), 4, 40)) for L in lengths]
    CRYSTALS[name] = dict(build=lambda a=atoms: a.copy(), sc=sc, mesh=mesh, normals={m: list(nvec)})
    if len(set(atoms.get_chemical_symbols())) > 1:
        import warnings
        warnings.warn(f"{name}: no Born charges given; LO-TO splitting omitted (add 'nac' in materials.CRYSTALS if polar)")
    return name, m


# ---------------------------------------------------------------- valence classes
_GROUP = {**{e: "II" for e in ("Zn", "Cd", "Mg", "Hg", "Be")}, **{e: "III" for e in ("B", "Al", "Ga", "In")},
          **{e: "IV" for e in ("C", "Si", "Ge", "Sn")}, **{e: "V" for e in ("N", "P", "As", "Sb")},
          **{e: "VI" for e in ("O", "S", "Se", "Te")}}


@__import__('functools').lru_cache(maxsize=None)
def compound_class(name):
    """'IV', 'III-V', 'II-VI', ... for covalent/ionic semiconductors; None otherwise."""
    if name not in CRYSTALS or name in METALS or name in SOMMERFELD:
        return None
    try:
        el = sorted(set(CRYSTALS[name]["build"]().get_chemical_symbols()))
    except Exception:
        return None
    g = sorted({_GROUP.get(e) for e in el}, key=lambda x: (x is None, x))
    if None in g:
        return None
    return "-".join(sorted(g, key=lambda x: ["II", "III", "IV", "V", "VI"].index(x)))


def heterovalent(a, b):
    """Label like 'III-V/IV' if the two semiconductors have different valence classes
    and one of them is elemental group IV (polar-on-nonpolar); else ''."""
    ca, cb = compound_class(a), compound_class(b)
    if ca and cb and ca != cb and "IV" in (ca, cb):
        return f"{ca}/{cb}"
    return ""
