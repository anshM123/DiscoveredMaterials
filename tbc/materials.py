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

# Thin-film thermal conductivities used for series (bulk) resistances of
# interlayers. TiC: ~20 W/mK (bulk TiC 21-33 W/mK); PMMA: 0.20 W/mK (accepted
# bulk value, Sandell et al. 2020); SiO2: 1.38 W/mK (fused silica).
FILM_K = {"TiC": 20.0, "PMMA": 0.20, "SiO2": 1.38,
          # 2 nm adhesion films: values used by Sandell et al. (2020), Table 1
          "Ti": 8.2, "Ni": 52.0}


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
AMORPHOUS = {"PMMA": PMMA, "SiO2": SIO2}
