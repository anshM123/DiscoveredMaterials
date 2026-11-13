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
     lattice-impedance spring factor (polymer contacts: physisorbed on noble metals → computed stiffness;
     chemisorbed on reactive metals and oxides → bonded, classified by oxide formation enthalpy);
   * an anharmonic (inelastic) channel;
   * the temperature convention of a measured conductance: TDTR temperature drops become local-equilibrium
     drops at high T (a fraction f(T) = 1/(1+(T_L/T)⁸) of the Zeng–Chen correction; f(300 K) ≈ 0);
   * the metal's electron–phonon resistance (Majumdar–Reddy, with DFT lattice conductivities; measured g,
     or Allen's formula with the MACE phonon DOS for metals without a measurement);
   * for metal/metal interfaces, the electron diffuse-mismatch channel in parallel;
   * series layers (thin films, adhesion layers coupled to the top metal's electrons).
3. **Accurate tier (~1 CPU-hour)** — an atomistic Green's function with MACE force constants of the actual
   relaxed coherent interface. It reproduces published DFT-NEGF for GaN/AlN (281 vs ~300 MW m⁻² K⁻¹). Being
   harmonic and defect-free, it differs from measured epitaxial wafers by 35% on average (AlN/SiC +4%,
   AlN/GaN −43%, GaN/SiC −52%, TiN/MgO −44%).
4. **Six global constants, frozen in `params.json` (none per specimen):**
   * α and b of the inelastic channel, and T_L, the onset of local equilibrium in a measured temperature
     drop, all fitted to one dataset (Al/Si, Li et al. 2022);
   * the specular efficiency η = 0.74, and 0.44 for heterovalent III-V/group-IV pairs, fitted to the 14
     MACE-AGF references (leave-one-out 14% / 29%);
   * the glass diffuson weight w, fitted to Li's Al/native-SiO₂/Si curve.

## Implementation

`tbc.py <A/hkl> <B/hkl> [--agf] [--polycrystalline] [--local-equilibrium] [--no-eph] [--T a:b:c]` runs on
any solid–solid interface, and `tbc.py --stack "Al/111 | SiO2:1.5 | Si/100"` on multilayers:

* any of the 30 library materials, or any crystal given as a structure file (`file.cif:hkl`);
* amorphous solids (glasses with a MACE melt-quench from `glass.py`; otherwise from density and sound velocities);
* metal sides are detected automatically (electron–phonon resistance; electron channel for metal/metal pairs).

It prints G(T), R_K(T) = 1/G, the phonon conductance G_pp (compare with anharmonic NEMD; `--local-equilibrium`
for NEMD's temperature convention), its harmonic elastic part G_el (compare with harmonic AGF/NEGF), Debye AMM
and DMM, and an expected-accuracy statement for the interface class. Runtime on a 4-core laptop CPU:

* a new crystal's phonons take 1–5 min;
* a G(T) curve takes seconds;
* an AGF reference takes ~1 h.

There are 24 physics-invariant tests and a smoke test covering every interface class.

## Results (from `make_report.py`; details and figures in README)

| test | this framework | strongest baseline |
|---|---|---|
| Li 2022 Al/Si G(T), 82–705 K, specimen 1 (fit set) | 3.6%; 300 K 395 vs 381; slope >400 K 0.31 vs 0.43 | 18% (full-dispersion DMM + same e-ph) |
| Li 2022 Al/Si G(T), specimen 2 (not fitted; 19% below specimen 1 at 300 K) | 25% | 12% (full-dispersion DMM + same e-ph) |
| Li 2022 Al/GaN G(T) (not fitted; guard set in model selection) | 5.7%; 300 K 436 vs 425; slope 0.28 vs 0.39 | 13% (full-dispersion DMM + same e-ph) |
| Cheng 2020 epitaxial Al/sapphire G(T) (not fitted; guard set) | 5.8%; 300 K 346 vs 345 | 4.6% (Debye AMM + same e-ph) |
| Sandell 2020 Au/[Ti,Ni]/PMMA/Si (not fitted; used to choose the polymer-contact treatment) | 13.0% (Au +2%, Ti −10%, Ni −26%); 34% if PMMA's Debye cutoff counts only C and O | 39% (Debye DMM + same e-ph/series) |
| Chang 2019 TiC enhancement (not fitted; measured 1.48; bulk-terminated diamond and phonon-only TiC assumed) | 1.44 (7.2 for H-terminated diamond; higher if TiC electrons carry heat) | 1.10 (Debye DMM) / 1.84 (full-dispersion DMM) |
| 14 coherent pairs, G_el vs MACE-AGF at 300 K (η calibration set; two-class η) | 15%; leave-one-out 14% (10 pairs) / 29% (4 III-V/IV pairs) | 66% on the original five (full-dispersion DMM) |
| 9 first-round blind pairs vs MACE-AGF (committed beforehand with the single η = 0.63 of the time; since then part of the η calibration) | 31%: the seven non-III-V/IV pairs 18%; **heterovalent III-V/IV failed** (GaP/Si(100) +54%, GaAs/Ge(100) +104%) | 46% (full-dispersion DMM) |
| Round-6 fresh blind pairs for the two-class η (committed before the AGF runs) | see README section B | |
| 6 measured epitaxial semiconductor and nitride/oxide interfaces (not fitted) | 21% (harmonic elastic part alone 41%) | 34% (full-dispersion DMM); 35% (harmonic AGF, 4 pairs); 236% (Debye DMM) |
| Li 2022 Al/native-SiO₂/Si G(T) (fit set of the glass diffuson weight; 8 points, 1 parameter) | 1.8% | 85% (Debye DMM + same e-ph) |

## Known limitations (stated, not hidden)

* The >400 K rise of Li's Al/Si (specimen 1) and Al/GaN is ~70% reproduced (0.31 vs 0.43 and 0.28 vs 0.39
  MW m⁻² K⁻²) with room-temperature values unchanged. Matching it fully is possible (0.41) but costs Al/GaN and
  Chang, so it is not taken. Li's second Al/Si specimen rises only 0.10 (model 0.26) and Cheng's Al/sapphire
  is flat to 479 K (model 0.30): the measured slope is specimen-dependent.
* Glasses: SiO₂ uses a MACE melt-quenched glass (propagons) plus one diffuson coupling weight, fitted to a
  single curve and tested only on Sandell. Amorphous Al₂O₃ has its own MACE glass (weight transferred from SiO₂, not
  validated). Other glasses need their own melt-quench (`glass.py`), or fall back to an upper bound.
* Molecule–metal chemistry (chemisorption, the Ni/PMMA case) is outside MACE-OMAT-0's domain, so such
  contacts are taken as bonded. The contact-classification rule was chosen after seeing the Sandell
  residuals (no parameter was fitted). Polymer stacks are the least accurate class (~10–25%).
* Harmonic atomistic references (AGF/NEGF) and anharmonic NEMD differ by ~2× for GaN/AlN. The tool reports
  both the harmonic (G_el, AGF) and the anharmonic (G_pp) estimate. The inelastic channel that separates them is
  calibrated on Al/Si only and adds 20–110% for nearly matched semiconductor pairs; it is checked against one
  anharmonic reference only (MLIP-NEMD AlN/GaN 600–700 vs 501 in the same temperature convention).
* Metal/metal interfaces use the electron diffuse-mismatch model (not validated here; ~2× high for Al/Cu).
  Metallic ceramics (TiN, TiC) are treated as phonon-only solids, so metal/TiN contacts are a lower bound.
* Al/GaN, Al/sapphire and the AGF pairs were never fitted, but served as no-regression guards in model
  selection (and the inelastic channel's Θ choice was made with Chang in view): they are validation sets.
  Pre-registered AGF predictions are the strictly blind tests. In the first round (single η = 0.63) the
  non-III-V/IV pairs were mostly low and the two heterovalent III-V/IV (100) pairs failed (+54%, +104%): MACE
  relaxes their electron-deficient interfacial bonds to soft, stretched ones. That motivated the two-class η
  (H14), whose calibration therefore includes those pairs; round 6 tests it on five fresh pairs. The CLI still
  recommends the AGF tier for III-V/IV interfaces.
* Below ~100 K the fast tier underestimates harmonic AGF for stiff nitrides by 25–55% (spectral undersampling).
* MACE-OMAT-0 is under an academic licence. A commercial universal potential is a drop-in replacement.
