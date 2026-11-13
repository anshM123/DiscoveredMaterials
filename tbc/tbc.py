#!/usr/bin/env python
"""Predict thermal boundary conductance G(T) and resistance R_K(T)=1/G(T) of a
solid-solid interface with the frozen model, alongside the AMM and DMM.

Examples
--------
  python tbc.py Al/111 Si/111 --metal Al
  python tbc.py AlN/0001 GaN/0001 --T 100:800:50
  python tbc.py Au PMMA --metal Au --Kn 1.2e19
  python tbc.py my_crystal.cif:001 Si/100          # any crystal from a structure file
  python tbc.py --stack "Al/111 | SiO2:1.5 | Si/100" --metal Al   # multilayer stack

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
from predict import load_params, G_ep, bonding_for, is_metal, G_edmm
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


def stack_main(argv):
    """tbc.py --stack "Al/111 | SiO2:1.5 | Si/100" [--metal Al] [--T ...]
    Layers separated by '|'; a finite film is NAME[/hkl]:thickness_nm. The first
    and last layers are semi-infinite. Interfaces in series, thin-film
    resistances (glassy films with Cahill-Watson-Pohl T-dependence), electron-
    phonon resistance of a metallic top layer (detected automatically or
    --metal). A thin metal film under the top metal (adhesion layer) also takes
    heat from the top metal's electrons through its own coupling g*d, in
    parallel; metal/metal interfaces add the electron channel (electron DMM).
    This is the treatment validated on the Sandell stacks (predict.sandell_stack)."""
    from predict import film_k, G_edmm
    ap = argparse.ArgumentParser()
    ap.add_argument("--stack", required=True)
    ap.add_argument("--metal", default=None)
    ap.add_argument("--no-eph", action="store_true")
    ap.add_argument("--T", default="100:700:100")
    ap.add_argument("--polycrystalline", action="store_true")
    a = ap.parse_args(argv)
    layers = []
    for tok in a.stack.split("|"):
        tok = tok.strip()
        name, _, thick = tok.partition(":")
        nm, sa = parse(name)
        ensure(nm)
        layers.append((nm, sa, float(thick) * 1e-9 if thick else None))
    P = load_params()
    if ":" in a.T:
        t0, t1, dt = map(float, a.T.split(":")); Ts = np.arange(t0, t1 + 1e-9, dt)
    else:
        Ts = [float(a.T)]
    metal = None if a.no_eph else (a.metal or (layers[0][0] if is_metal(layers[0][0]) else None))
    adh = layers[1] if (metal and len(layers) > 2 and layers[1][2] and is_metal(layers[1][0])) else None
    ifs = []
    for (n1, s1_, _), (n2, s2_, _) in zip(layers[:-1], layers[1:]):
        bond = bonding_for(n1, n2, s1_, s2_)
        ifs.append((Interface(Side(n1, s1_), Side(n2, s2_), K=bond[0] if bond else None,
                              Kref=bond[1] if bond else None, epitaxial=not a.polycrystalline),
                    (n1, n2) if (is_metal(n1) and is_metal(n2)) else None))
    print("stack:", " | ".join(f"{n}{'/' + s if s else ''}{'(%.1f nm)' % (t * 1e9) if t else ''}" for n, s, t in layers),
          f"  metal: {metal or '-'}" + (f" (+ adhesion film {adh[0]} coupled in parallel)" if adh else ""))
    print("   T(K)   G_total(MW/m2K)  R_total(m2K/GW)   per-interface G (MW/m2K)")
    for T in Ts:
        R = 0.0
        if metal:
            gep = G_ep(metal, T) + (G_ep(adh[0], T, d=adh[2]) if adh else 0.0)
            R += 1 / gep
        gs = []
        for I, mm in ifs:
            g = I.G_pp(T, P["alpha"], P["b"])
            if mm:
                g += G_edmm(*mm, T)
            gs.append(g)
        R += sum(1 / g for g in gs)
        for n, _, t in layers[1:-1]:
            if t:
                R += t / film_k(n, T)
        print(f"{T:7.0f} {1e-6 / R:14.1f} {1e9 * R:15.3f}    " + "  ".join(f"{g / 1e6:.0f}" for g in gs))


def main():
    import sys as _sys
    if "--stack" in _sys.argv:
        return stack_main(_sys.argv[1:])
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--T", default="100:700:100", help="start:stop:step in K, or a single T")
    ap.add_argument("--metal", default=None, help="name of the metallic side (adds electron-phonon resistance); "
                    "elemental metals are detected automatically")
    ap.add_argument("--no-eph", action="store_true", help="phonon conductance only (e.g. to compare with "
                    "atomistic MD/AGF, which have no electrons)")
    ap.add_argument("--Kn", type=float, default=None, help="interfacial normal stiffness, N/m^3")
    ap.add_argument("--polycrystalline", action="store_true",
                    help="non-epitaxial (sputtered/polycrystalline) interface: fully diffuse scattering")
    ap.add_argument("--agf", action="store_true",
                    help="accurate tier for coherent crystalline pairs: atomistic Green's function with MACE "
                         "force constants (~1 CPU-hour; harmonic, perfect interface). Cached in data/agf.")
    ap.add_argument("--local-equilibrium", action="store_true",
                    help="also print the phonon conductance in the local-equilibrium-temperature convention "
                         "(Zeng & Chen correction), the convention of NEMD temperature drops")
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
    agf_res = None
    if a.agf:
        import agf as agfmod
        lab = f"{na}{sa}_{nb}{sb}"
        path = os.path.join(agfmod.CACHE, f"{lab}.json")
        if not os.path.exists(path):
            print(f"[running AGF for {lab}; this takes ~1 CPU-hour]")
            agfmod.run(na, sa, nb, sb, nk=8, nw=160, inplane_min=9.0)
        agf_res = json.load(open(path))
    metals = [] if a.no_eph else ([a.metal] if a.metal else [n for n in (na, nb) if is_metal(n)])
    if a.metal and not a.no_eph and is_metal(na) and is_metal(nb):
        metals = [na, nb]
    rows = []
    for T in Ts:
        Gpp = I.G_pp(T, P["alpha"], P["b"])
        Ge = 0.0
        if len(metals) == 2:
            # metal/metal: electrons cross directly (electron DMM), in parallel with
            # the phonon path e-ph(1) -> interface phonons -> e-ph(2)
            Ge = G_edmm(metals[0], metals[1], T)
            G = Ge + 1 / (1 / G_ep(metals[0], T) + 1 / Gpp + 1 / G_ep(metals[1], T))
        elif metals:
            G = 1 / (1 / Gpp + 1 / G_ep(metals[0], T))
        else:
            G = Gpp
        rows.append(dict(T=T, G=G / 1e6, R_K=1e9 / G, G_pp=Gpp / 1e6, G_electron=Ge / 1e6,
                         G_pp_local_eq=I.G_pp(T, P["alpha"], P["b"], noneq=True) / 1e6,
                         G_el_local_eq=I.G_pp(T, 0.0, 1.0, noneq=True) / 1e6,
                         G_elastic=I.G_elastic(T) / 1e6, P_inelastic=I.P_inel(T, P["alpha"], P["b"]),
                         AMM=G_AMM_debye(s1, s2, T) / 1e6, DMM=G_DMM_debye(s1, s2, T) / 1e6))
        if agf_res is not None:
            g_agf = float(np.interp(T, agf_res["T"], agf_res["G"]))
            rows[-1]["G_AGF_harmonic"] = g_agf / 1e6
            if len(metals) == 1:
                rows[-1]["G_AGF_with_eph"] = 1 / (1 / g_agf + 1 / G_ep(metals[0], T)) / 1e6
    if a.json:
        fin = lambda v: (None if isinstance(v, float) and not np.isfinite(v) else v)
        rows_j = [{k: fin(float(v)) for k, v in r.items()} for r in rows]
        print(json.dumps(dict(interface=f"{a.a}|{a.b}", bonding=bond, metals=metals,
                              units="G in MW m^-2 K^-1, R_K in m^2 K GW^-1; null = infinite (no interface)",
                              rows=rows_j), indent=1)); return
    print(f"interface {a.a} | {a.b}   bonding: {'ideal' if bond is None else 'K/Kref = %.2f' % (bond[0][0] / bond[1][0])}"
          f"   metal: {'/'.join(metals) or '-'}   coherence length: {'diffuse' if I.Lc is None else '%.1f A' % I.Lc}")
    # expected accuracy of this interface class, from the validation in README
    if not (s1.crystal and s2.crystal):
        import model as _model
        amo = [x.name for x in (s1, s2) if not x.crystal]
        if any(n in _model.GLASS_FLUX for n in amo):
            acc = ("glass with a MACE glass model (propagons + fitted diffuson weight): Al/native-SiO2/Si G(T) "
                   "within ~2% (fit set); metal/PMMA/SiO2 stacks ~12% mean (Ni -25%)")
        elif "PMMA" in amo:
            acc = "polymer (Debye solid): metal/PMMA/SiO2/Si stacks ~12% mean error (Au +2%, Ti -9%, Ni -25%)"
        else:
            acc = ("amorphous side without a MACE glass model (Debye glass): treat as an UPPER BOUND "
                   "(Debye SiO2 over-predicted Al/native-SiO2/Si by ~1.5x); build one with glass.py")
    elif len(metals) == 2:
        acc = ("metal/metal: electrons cross directly (electron DMM, column G_e); EDMM over-predicts Al/Cu ~2x "
               "(Gundrum et al. 2005, ~4 GW measured). Not validated here: factor-2 estimate")
    elif I.Lc is not None and I.p_inf > 0.5:
        acc = ("coherent crystalline pair: G_el within ~19% of harmonic AGF (leave-one-out, 5 pairs; 3 blind tests "
               "-13/-12/-7%); measured epitaxial interfaces ~25% mean. For atomistic accuracy run --agf")
    elif metals or any(is_metal(n) for n in (na, nb)):
        acc = ("metal on crystal (incommensurate/diffuse): validated for Al on Si, GaN, sapphire (5-7%; a second "
               "Al/Si specimen 24%, specimen-to-specimen spread ~20%); other metals not validated")
    else:
        acc = "incommensurate or polycrystalline crystalline pair: diffuse limit; validated only indirectly (~25%)"
    print(f"expected accuracy: {acc}")
    print("columns: G = total (phonons incl. inelastic channel, electron-phonon resistance of metal sides, electron "
          "channel for metal/metal); G_pp = phonon interface conductance (compare with anharmonic NEMD); G_el = its "
          "harmonic elastic part (compare with harmonic AGF/NEGF)")
    hdr = "   T(K)   G(MW/m2K)  R_K(m2K/GW)  G_pp   G_el  P_inel   AMM    DMM"
    print(hdr + ("     G_e" if len(metals) == 2 else "") + ("   G_pp(LE)  G_el(LE)" if a.local_equilibrium else ""))
    for r in rows:
        print("{T:7.0f} {G:10.1f} {R_K:11.3f} {G_pp:7.0f} {G_elastic:6.0f} {P_inelastic:6.2f} {AMM:6.0f} {DMM:6.0f}".format(**r)
              + (" {:7.0f}".format(r["G_electron"]) if len(metals) == 2 else "")
              + (" {:10.0f} {:9.0f}".format(r["G_pp_local_eq"], r["G_el_local_eq"]) if a.local_equilibrium else "")
              + ("   AGF(harmonic) {:.0f}".format(r["G_AGF_harmonic"]) if "G_AGF_harmonic" in r else ""))


if __name__ == "__main__":
    main()
