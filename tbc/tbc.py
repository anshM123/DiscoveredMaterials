#!/usr/bin/env python
"""Predict thermal boundary conductance G(T) and resistance R_K(T)=1/G(T) of a
solid-solid interface with the frozen model, alongside the AMM and DMM.

Examples
--------
  python tbc.py Al/111 Si/111 --metal Al
  python tbc.py AlN/0001 GaN/0001 --T 100:800:50
  python tbc.py Au PMMA --metal Au --Kn 1.2e19
  python tbc.py my_crystal.cif:001 Si/100          # any crystal from a structure file

A material must be defined in materials.py (CRYSTALS: a structure guess; or
AMORPHOUS: density, sound velocities). Phonons of a new crystal are computed
with MACE-OMAT-0 on first use and cached (a few minutes on a laptop CPU).

Bonding: by default the interface is treated as chemically (bulk-like)
bonded. For weakly bonded interfaces pass the interfacial normal stiffness
--Kn (N/m^3), e.g. from bonding.py (MACE+D3 rigid-shift scan) or from DFT.
"""
import argparse, json, os
import numpy as np
from phonons import compute, CACHE
from materials import CRYSTALS
from model import Side, Interface, G_AMM_debye, G_DMM_debye
from predict import load_params, G_ep, bonding_for
from bonding import Kref


def parse(spec):
    """'Name/normal' (library material) or 'path/to/file.cif:hkl' (any crystal)."""
    if ":" in spec and os.path.exists(spec.rsplit(":", 1)[0]):
        from materials import register_structure
        path, hkl = spec.rsplit(":", 1)
        return register_structure(path, hkl)
    name, _, normal = spec.partition("/")
    from materials import AMORPHOUS
    if name not in CRYSTALS and name not in AMORPHOUS:
        raise SystemExit(f"error: unknown material '{name}'. Library: {sorted(CRYSTALS) + sorted(AMORPHOUS)}.\n"
                         f"Pass a structure file instead, e.g. path/to/{name}.cif:111, or add it to materials.py.")
    if normal and name in CRYSTALS and normal not in CRYSTALS[name]["normals"]:
        raise SystemExit(f"error: normal '{normal}' not available for {name}; options: {list(CRYSTALS[name]['normals'])}")
    return name, (normal or None)


def ensure(name):
    if name in CRYSTALS and not os.path.exists(os.path.join(CACHE, f"{name}.npz")):
        print(f"[computing phonons of {name} with MACE-OMAT-0 ...]")
        compute(name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--T", default="100:700:100", help="start:stop:step in K, or a single T")
    ap.add_argument("--metal", default=None, help="name of the metallic side (adds electron-phonon resistance)")
    ap.add_argument("--Kn", type=float, default=None, help="interfacial normal stiffness, N/m^3")
    ap.add_argument("--polycrystalline", action="store_true",
                    help="non-epitaxial (sputtered/polycrystalline) interface: fully diffuse scattering")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    (na, sa), (nb, sb) = parse(a.a), parse(a.b)
    ensure(na); ensure(nb)
    P = load_params()
    if a.Kn is not None:
        kr = Kref(na, nb, sa, sb); r = a.Kn / kr["K_n"]
        bond = ((r * kr["K_n"], r * kr["K_t"]), (kr["K_n"], kr["K_t"]))
    else:
        bond = bonding_for(na, nb, sa, sb)
    s1, s2 = Side(na, sa), Side(nb, sb)
    I = Interface(s1, s2, K=bond[0] if bond else None, Kref=bond[1] if bond else None,
                  epitaxial=not a.polycrystalline)
    if ":" in a.T:
        t0, t1, dt = map(float, a.T.split(":")); Ts = np.arange(t0, t1 + 1e-9, dt)
    else:
        Ts = [float(a.T)]
    rows = []
    for T in Ts:
        Gpp = I.G_pp(T, P["alpha"], P["b"])
        G = 1 / (1 / Gpp + 1 / G_ep(a.metal, T)) if a.metal else Gpp
        rows.append(dict(T=T, G=G / 1e6, R_K=1e9 / G, G_pp=Gpp / 1e6,
                         G_elastic=I.G_elastic(T) / 1e6, P_inelastic=I.P_inel(T, P["alpha"], P["b"]),
                         AMM=G_AMM_debye(s1, s2, T) / 1e6, DMM=G_DMM_debye(s1, s2, T) / 1e6))
    if a.json:
        print(json.dumps(dict(interface=f"{a.a}|{a.b}", bonding=bond, units="G in MW m^-2 K^-1, R_K in m^2 K GW^-1",
                              rows=rows), indent=1)); return
    print(f"interface {a.a} | {a.b}   bonding: {'ideal' if bond is None else 'K/Kref = %.2f' % (bond[0][0] / bond[1][0])}"
          f"   metal: {a.metal or '-'}   coherence length: {'diffuse' if I.Lc is None else '%.1f A' % I.Lc}")
    print("   T(K)   G(MW/m2K)  R_K(m2K/GW)  G_pp   G_el  P_inel   AMM    DMM")
    for r in rows:
        print("{T:7.0f} {G:10.1f} {R_K:11.3f} {G_pp:7.0f} {G_elastic:6.0f} {P_inelastic:6.2f} {AMM:6.0f} {DMM:6.0f}".format(**r))


if __name__ == "__main__":
    main()
