### A. Experiments (fast tier)

Mean absolute error against measurement. AMM/DMM are the Debye models of Little (1959) and Swartz & Pohl (1989), shown raw and with the same electron–phonon resistance and series layers as this model (`+e-ph`). *Guard set: never fitted, but used as a no-regression guard when accepting or rejecting model changes (EXPERIMENTS.md), so it is a validation set rather than a blind test set. The last column is the strongest baseline: full-dispersion DMM on the same MACE phonon spectra plus the same e-ph resistance (no specularity, bonding or inelastic channel). Only Al/Si sample 1 is used for fitting (2 global parameters).

| dataset (role) | this model | AMM | DMM | AMM+e-ph | DMM+e-ph | full-disp. DMM+e-ph |
|---|---|---|---|---|---|---|
| Al/Si sample1, 80–700 K (calibration) | 4.6% | 172.8% | 105.7% | 57.8% | 32.7% | 18.2% |
| Al/Si sample2, 80–700 K (second specimen, not fitted) | 24.3% | 238.9% | 155.8% | 100.2% | 67.1% | 12.4% |
| Al/GaN, 80–700 K (not fitted; guard set*) | 6.9% | 180.9% | 128.3% | 54.0% | 36.8% | 12.9% |
| Al/sapphire (Cheng 2020), 80–700 K (not fitted; guard set*) | 6.2% | 32.5% | 141.8% | 4.6% | 54.1% | 14.1% |
| Sandell 2020, 3 metal/PMMA/Si stacks (not fitted; used to choose the contact treatment, see EXPERIMENTS H8–H9) | 13.0% | 51.2% | 79.4% | 47.0% | 39.2% | – |
| Chang 2019 TiC enhancement ratio, measured 1.48 (not fitted; bulk-terminated diamond assumed) | 1.36 | 2.19 | 1.10 | 2.00 | 1.08 | 1.84 |

The two Al/Si specimens of Li et al. differ by 23% at 300 K (381 vs 310 MW m⁻² K⁻¹): specimen-to-specimen variation of nominally identical interfaces is itself ~20–25%, which bounds the accuracy any structure-only model can claim for real samples.

Sandell stacks (MW m⁻² K⁻¹, measured / model): Au/PMMA/Si 59 / 60, Au/Ti/PMMA/Si 115 / 103, Au/Ni/PMMA/Si 139 / 102

### B. Coherent crystalline interfaces: fast tier vs MACE-AGF (300 K, harmonic)

The AGF tier (`agf.py`, `tbc.py --agf`) is the reference for perfect coherent interfaces. Independent check: AlN/GaN 281 vs ~300 MW m⁻² K⁻¹ from DFT-NEGF (Polanco & Lindsay, PRB 99, 075202, 2019). Harmonic AGF describes a perfect interface without anharmonic channels; against measured wafers it is −43% (AlN/GaN), +4% (AlN/SiC) and −52% (GaN/SiC) (section C). **Which fast-tier output to compare with which reference:** harmonic AGF/NEGF ↔ `G_el`; anharmonic NEMD or experiment ↔ `G` (adds the inelastic channel calibrated on Al/Si). The table gives both; `G` is expected to exceed harmonic AGF. `G_el (CLI)` is the CLI's automatic coincidence-matching path, `G_el` the coherent setting used to calibrate η.

| interface | strain in AGF | **AGF** | fast tier G_el | G_el (CLI) | fast tier G (anharmonic) | full-dispersion DMM | DMM (Debye) | AMM (Debye) |
|---|---|---|---|---|---|---|---|---|
| AlN(0001)/GaN(0001) | -2.9% | **281** | 243 (-14%) | 246 (-12%) | 378 | 445 (+59%) | 1376 (+390%) | 1529 (+444%) |
| AlN(0001)/SiC4H(0001) | +1.1% | **487** | 531 (+9%) | 532 (+9%) | 657 | 755 (+55%) | 1849 (+280%) | 2844 (+484%) |
| GaAs(100)/AlAs(100) | +0.7% | **170** | 133 (-21%) | 133 (-21%) | 277 | 171 (+1%) | 594 (+250%) | 928 (+446%) |
| GaN(0001)/SiC4H(0001) | +4.1% | **147** | 184 (+25%) | 188 (+28%) | 311 | 374 (+154%) | 1291 (+776%) | 1084 (+636%) |
| Si(100)/Ge(100) | -5.2% | **126** | 135 (+8%) | 137 (+9%) | 254 | 204 (+62%) | 547 (+335%) | 403 (+221%) |
| **mean abs. deviation from AGF (calibration set of η)** | | | 15% | 16% | – | 66% | 406% | 446% |

**Blind tests** (fast-tier value committed before the AGF reference was computed; η frozen):

| interface | **AGF** | pre-registered fast tier G_el | G_el (CLI, now) | full-dispersion DMM | DMM (Debye) | AMM (Debye) |
|---|---|---|---|---|---|---|
| GaN(0001)/ZnO(0001) | **207** | 181 (-13%) | 182 (-12%) | 307 (+48%) | 1099 (+430%) | 1363 (+557%) |
| GaAs(110)/AlAs(110) | **154** | 135 (-12%) | 135 (-12%) | 166 (+8%) | 594 (+286%) | 928 (+504%) |
| TiN(100)/MgO(100) | **394** | 366 (-7%) | 366 (-7%) | 470 (+19%) | 2167 (+450%) | 3971 (+907%) |
| GaP(100)/Si(100) | **105** | 163 (+56%) | 163 (+56%) | 249 (+138%) | 760 (+627%) | 976 (+833%) |

### C. Measured conductance of epitaxial semiconductor and nitride/oxide interfaces (300 K, not fitted)

Multi-frequency TDTR on GaN-on-SiC, GaN-on-AlN and AlN-on-SiC wafers (par.nsf.gov 10498482). * Epitaxial TiN on MgO(001), MgO(111) and sapphire: Costescu, Wall & Cahill, PRB 67, 054302 (2003) found the same conductance on all three; ~700 MW m⁻² K⁻¹ at room temperature as cited by Wilson et al., PRB 91, 115414 (2015) (approximate; TiN electron-phonon term omitted: no g data). Real interfaces contain defects and anharmonic transport; harmonic AGF is a perfect-interface reference.

| interface | measured | fast tier (total) | fast tier (elastic) | AGF (harmonic) | DMM (Debye) | AMM (Debye) |
|---|---|---|---|---|---|---|
| AlN(0001)/GaN(0001) | 490 | 378 | 246 | 281 | 1376 | 1529 |
| AlN(0001)/4H-SiC(0001) | 470 | 657 | 532 | 487 | 1849 | 2844 |
| GaN(0001)/4H-SiC(0001) | 305 | 311 | 188 | 147 | 1291 | 1084 |
| TiN(001)/MgO(001) * | 700 | 512 | 366 | 394 | 2167 | 3971 |
| TiN(111)/MgO(111) * | 700 | 496 | 338 | – | 2167 | 3971 |
| TiN(111)/sapphire(0001) * | 700 | 482 | 372 | – | 2088 | 3620 |
| **mean abs. error** | | 25% | | 35% | 236% | 387% |
| mean abs. error, the 4 interfaces with an AGF reference | | 23% | | 35% | 252% | 360% |

### D. Metal on amorphous native oxide: Al / SiO₂(native) / Si, 80–350 K (calibration set of the glass diffuson weight)

Li et al. 2022, Fig. 1a (orange markers). The native-oxide thickness is not reported; 1.0–2.0 nm is typical (Sandell et al. measured 2.05 nm on their wafer). The film conductivity follows the Cahill–Watson–Pohl T-dependence anchored at 1.38 W m⁻¹ K⁻¹ (300 K).

| model | mean abs. error | G(80 K) | G(300 K) |
|---|---|---|---|
| measured | – | 52 | 138 |
| this model, oxide 1.5 nm (fit of w) | 2% | 53 | 140 |
| this model, oxide 1.0 nm | 6% | 56 | 147 |
| this model, oxide 2.0 nm | 6% | 50 | 133 |
| DMM (Debye) + e-ph, oxide 1.5 nm | 85% | 91 | 266 |
| AMM (Debye) + e-ph, oxide 1.5 nm | 128% | 118 | 319 |

The glass is described by the propagating flux of a MACE melt-quenched a-SiO₂ plus a fraction w of its diffuson channels; w = 0.27 is fitted to this curve (oxide 1.5 nm; 1.0–2.0 nm gives w = 0.24–0.32). One parameter reproduces the whole 80–350 K temperature dependence. The Sandell stacks (PMMA on native SiO₂/Si) were not used in the fit of w (they were used to choose the polymer-contact treatment, EXPERIMENTS H8–H9).

![parity](results/parity.png)

*Parity plot of every validated case: experiments (Li, Cheng, Sandell; all temperatures), measured epitaxial interfaces (section C) and MACE-AGF references (section B), for this model and for the Debye AMM/DMM.*
