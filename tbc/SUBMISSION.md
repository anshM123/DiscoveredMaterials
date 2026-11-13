# Submission — Grand Challenge 05: Predicting Heat Transport across Solid Interfaces

**To:** the submission address on the Grand Challenges page.
**Code, data and full documentation:** this repository, directory `tbc/` (README.md). Everything
reproduces with `pip install -r requirements.txt && python make_report.py`.

## Method

A two-tier, structure-to-conductance framework. Nothing is fitted per specimen.

1. **Phonons from a universal machine-learned potential.** MACE-OMAT-0 force constants and phonopy give
   full-dispersion, interface-normal phonon fluxes, including the in-plane-wavevector-resolved channel
   counts. Polar crystals include LO–TO splitting.
2. **Fast tier (seconds per G(T) curve)** — a Landauer model with:
   * momentum-conserving specular transmission for coherent interfaces, plus diffuse (full-dispersion DMM)
     transmission, mixed by a specularity computed from the automatic coincidence-lattice match of the
     two interface planes;
   * interfacial bond stiffness from MACE+D3 rigid-shift scans of the actual interface, entering through a
     lattice-impedance spring factor;
   * an anharmonic (inelastic) channel;
   * the metal's electron–phonon resistance (Majumdar–Reddy, with DFT lattice conductivities);
   * series layers.
3. **Accurate tier (~1 CPU-hour)** — an atomistic Green's function with MACE force constants of the actual
   relaxed coherent interface. It reproduces published DFT-NEGF for GaN/AlN (281 vs ~300 MW m⁻² K⁻¹) and the
   measured AlN/4H-SiC (487 vs ~470).
4. **Three global parameters, frozen in `params.json`:**
   * α and b of the inelastic channel, fitted to one dataset (Al/Si, Li et al. 2022);
   * the specular efficiency η, fitted to five MACE-AGF references.

## Implementation

`tbc.py <A/hkl> <B/hkl> [--metal M] [--agf] [--polycrystalline] [--T a:b:c]` runs on any solid–solid
interface:

* any of the 30 library materials, or any crystal given as a structure file (`file.cif:hkl`);
* amorphous solids from their density and sound velocities.

It prints G(T), R_K(T) = 1/G, the elastic (harmonic) and total conductance, Debye AMM and DMM, and an
expected-accuracy statement for the interface class. Runtime on a 4-core laptop CPU:

* a new crystal's phonons take 1–5 min;
* a G(T) curve takes seconds;
* an AGF reference takes ~1 h.

There are 19 physics-invariant tests and a smoke test covering every interface class.

## Results (from `make_report.py`; details and figures in README)

| test | this framework | best of Debye AMM/DMM |
|---|---|---|
| Li 2022 Al/Si G(T), 80–700 K (fit set) | 4.6% | 106% |
| Li 2022 Al/GaN G(T) (held out) | 6.9% | 128% |
| Cheng 2020 epitaxial Al/sapphire G(T) (held out) | 6.2% | 77% |
| Sandell 2020 Au/[Ti,Ni]/PMMA/Si (held out) | 19.6% | 22% (AMM + this e-ph/series treatment: 10% — better than this model) |
| Chang 2019 TiC enhancement (held out; measured 1.48) | 1.36 | 1.97 / 1.10 |
| 5 coherent pairs vs MACE-AGF (η fit set; leave-one-out) | 15% (LOO 19%) | 406% |
| 3 blind coherent pairs vs MACE-AGF (predictions committed beforehand) | 11% (−13, −12, −7%) | 389% |
| 6 measured epitaxial semiconductor and nitride/oxide interfaces (held out) | 25% | 236% |
| Li 2022 Al/native-SiO₂/Si G(T) (held out; oxide 1–2 nm) | 36–60% (**weak**) | 113% |

## Known limitations (stated, not hidden)

* The >400 K rise of Li's Al/Si and Al/GaN is reproduced in sign but about 4× too weak (0.10 vs 0.43 MW m⁻² K⁻²).
  Calibrating to the slope reproduces it but costs the room-temperature values (trade-off table in the
  README). Epitaxial Al/sapphire (Cheng) is flat, and no single global mechanism fits both groups' data.
* Glasses are Debye solids: a metal on native SiO₂ is over-predicted ~1.4–1.7×, so metal/amorphous-dielectric
  results are upper bounds. An atomistic glass model (MACE melt-quench + AGF) would be the fix.
* Molecule–metal chemistry, the Ni/PMMA case, is outside MACE-OMAT-0's domain. Polymer contacts are
  therefore the least accurate class (~20–40%).
* Harmonic atomistic references (AGF/NEGF) and anharmonic NEMD differ by ~2× for GaN/AlN. The tool reports
  both the harmonic (G_elastic, AGF) and the anharmonic (G) estimate.
* MACE-OMAT-0 is under an academic licence. A commercial universal potential is a drop-in replacement.
