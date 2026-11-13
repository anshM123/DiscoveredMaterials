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
   * the metal's electron–phonon resistance (Majumdar–Reddy, with DFT lattice conductivities; measured g,
     or Allen's formula with the MACE phonon DOS for metals without a measurement);
   * for metal/metal interfaces, the electron diffuse-mismatch channel in parallel;
   * series layers (thin films, adhesion layers coupled to the top metal's electrons).
3. **Accurate tier (~1 CPU-hour)** — an atomistic Green's function with MACE force constants of the actual
   relaxed coherent interface. It reproduces published DFT-NEGF for GaN/AlN (281 vs ~300 MW m⁻² K⁻¹). Being
   harmonic and defect-free, it differs from measured epitaxial wafers by 35% on average (AlN/SiC +4%,
   AlN/GaN −43%, GaN/SiC −52%).
4. **Four global parameters, frozen in `params.json`:**
   * α and b of the inelastic channel, fitted to one dataset (Al/Si, Li et al. 2022);
   * the specular efficiency η, fitted to five MACE-AGF references;
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

There are 23 physics-invariant tests and a smoke test covering every interface class.

## Results (from `make_report.py`; details and figures in README)

| test | this framework | strongest baseline |
|---|---|---|
| Li 2022 Al/Si G(T), 80–700 K, specimen 1 (fit set) | 4.6% | 18% (full-dispersion DMM + same e-ph) |
| Li 2022 Al/Si G(T), specimen 2 (not fitted; 23% below specimen 1 at 300 K) | 24% | 12% (full-dispersion DMM + same e-ph) |
| Li 2022 Al/GaN G(T) (not fitted; guard set in model selection) | 6.9% | 13% (full-dispersion DMM + same e-ph) |
| Cheng 2020 epitaxial Al/sapphire G(T) (not fitted; guard set) | 6.2% | 4.6% (Debye AMM + same e-ph) |
| Sandell 2020 Au/[Ti,Ni]/PMMA/Si (not fitted; used to choose the polymer-contact treatment) | 13.0% (Au +2%, Ti −10%, Ni −26%); 34% if PMMA's Debye cutoff counts only C and O | 39% (Debye DMM + same e-ph/series) |
| Chang 2019 TiC enhancement (not fitted; measured 1.48; bulk-terminated diamond and phonon-only TiC assumed) | 1.36 (1.9–3.2 if TiC electrons carry heat) | 1.10 (Debye DMM) / 1.84 (full-dispersion DMM) |
| 5 coherent pairs, G_el vs MACE-AGF at 300 K (η fit set; leave-one-out) | 15% (LOO 19%); 28% at 50 K | 66% (full-dispersion DMM) |
| 4 blind coherent pairs vs MACE-AGF (predictions committed beforehand) | 22% (−13, −12, −7, **+54%** for polar GaP/Si(100)) | 53% (full-dispersion DMM) |
| 6 measured epitaxial semiconductor and nitride/oxide interfaces (not fitted) | 25% | 35% (harmonic AGF); 236% (Debye DMM) |
| Li 2022 Al/native-SiO₂/Si G(T) (fit set of the glass diffuson weight; 8 points, 1 parameter) | 1.8% | 85% (Debye DMM + same e-ph) |

## Known limitations (stated, not hidden)

* The >400 K rise of Li's Al/Si (specimen 1) and Al/GaN is reproduced in sign but about 4× too weak (0.10 vs
  0.43 MW m⁻² K⁻²). Li's second Al/Si specimen rises only 0.10 (model 0.14) and Cheng's Al/sapphire is flat
  (model 0.25): the measured slope is specimen-dependent. Calibrating to the slope costs the room-temperature
  values (trade-off table in the README).
* Glasses: SiO₂ uses a MACE melt-quenched glass (propagons) plus one diffuson coupling weight, fitted to a
  single curve and tested only on Sandell. Other glasses need their own melt-quench, or fall back to an upper bound.
* Molecule–metal chemistry (chemisorption, the Ni/PMMA case) is outside MACE-OMAT-0's domain, so such
  contacts are taken as bonded. The contact-classification rule was chosen after seeing the Sandell
  residuals (no parameter was fitted). Polymer stacks are the least accurate class (~10–25%).
* Harmonic atomistic references (AGF/NEGF) and anharmonic NEMD differ by ~2× for GaN/AlN. The tool reports
  both the harmonic (G_el, AGF) and the anharmonic (G_pp) estimate. The inelastic channel that separates them is
  calibrated on Al/Si only and adds 20–110% for nearly matched semiconductor pairs; it is unverified against NEMD.
* Metal/metal interfaces use the electron diffuse-mismatch model (not validated here; ~2× high for Al/Cu).
  Metallic ceramics (TiN, TiC) are treated as phonon-only solids, so metal/TiN contacts are a lower bound.
* Al/GaN, Al/sapphire and the AGF pairs were never fitted, but served as no-regression guards in model
  selection (and the inelastic channel's Θ choice was made with Chang in view): they are validation sets.
  The four pre-registered AGF predictions are the strictly blind tests (one failed: GaP/Si(100) +54%).
* Below ~100 K the fast tier underestimates harmonic AGF for stiff nitrides by 25–55% (spectral undersampling).
* MACE-OMAT-0 is under an academic licence. A commercial universal potential is a drop-in replacement.
