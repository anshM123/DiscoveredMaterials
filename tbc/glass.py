"""Amorphous SiO2 from MACE-OMAT-0 melt-quench, in a cell laterally matched to
Al(111) 4x4 and Si(111) 3x3 (hexagonal, a = A_LAT), for atomistic Green's-function
references of metal/glass and glass/semiconductor interfaces.

Protocol: random SiO2 at the experimental density 2.20 g/cm^3 -> Langevin MD at
5000 K (4 ps) -> linear quench to 300 K (12 ps) -> 300 K (2 ps) -> relaxation
(float64). Cached in data/glass/.
"""
import os, sys, json, warnings
import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "glass")
A_LAT = 11.45            # Angstrom: Al(111) 4x4 = 11.41, Si(111) 3x3 = 11.51
RHO = 2.20               # g/cm^3


def make_glass(n_units=32, seed=0, T_melt=5000.0, t_melt=4.0, t_quench=12.0, t_anneal=2.0, dt=1.0):
    from ase import Atoms, units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.optimize import FIRE
    from bonding import calc32, calc
    rng = np.random.default_rng(seed)
    mass = n_units * 60.084 / 6.02214076e23          # g
    vol = mass / RHO * 1e24                            # A^3
    area = A_LAT ** 2 * np.sqrt(3) / 2
    c = vol / area
    cell = [[A_LAT, 0, 0], [-A_LAT / 2, A_LAT * np.sqrt(3) / 2, 0], [0, 0, c]]
    # random placement with a minimum distance
    pos, sym = [], []
    for s_, n in (("Si", n_units), ("O", 2 * n_units)):
        for _ in range(n):
            for _try in range(2000):
                f = rng.random(3)
                p = f @ np.array(cell)
                ok = True
                if pos:
                    a = Atoms(sym + [s_], positions=pos + [p], cell=cell, pbc=True)
                    d = a.get_distances(len(a) - 1, list(range(len(a) - 1)), mic=True)
                    ok = d.min() > 1.6
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


def analyse(at):
    """Coordination statistics: fraction of 4-fold Si and 2-fold O."""
    from ase.neighborlist import neighbor_list
    i, j = neighbor_list("ij", at, 2.0)
    sy = np.array(at.get_chemical_symbols())
    cn = np.bincount(i, minlength=len(at))
    si = sy == "Si"; o = sy == "O"
    return dict(Si4=float(np.mean(cn[si] == 4)), O2=float(np.mean(cn[o] == 2)),
                density=float(np.sum(at.get_masses()) / 6.02214076e23 / at.get_volume() * 1e24))


if __name__ == "__main__":
    from ase.io import write
    os.makedirs(OUT, exist_ok=True)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    g = make_glass(n_units=n)
    st = analyse(g)
    print(json.dumps(st))
    write(os.path.join(OUT, f"aSiO2_{n}.extxyz"), g)
    json.dump(st, open(os.path.join(OUT, f"aSiO2_{n}.json"), "w"))
