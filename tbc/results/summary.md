### A. Experiments (fast tier)

Mean absolute error against measurement. AMM/DMM are the Debye models of Little (1959) and Swartz & Pohl (1989), shown raw and with the same electron–phonon resistance and series layers as this model (`+e-ph`). Only Al/Si is used for fitting (2 global parameters).

| dataset (role) | this model | AMM | DMM | AMM+e-ph | DMM+e-ph |
|---|---|---|---|---|---|
| Al/Si sample1, 80–700 K (calibration) | 4.6% | 183.0% | 105.7% | 61.6% | 32.7% |
| Al/GaN, 80–700 K (**held out**) | 6.9% | 225.2% | 128.3% | 68.3% | 36.8% |
| Al/sapphire (Cheng 2020), 80–700 K (**held out**) | 6.2% | 76.7% | 141.8% | 27.0% | 54.1% |
| Sandell 2020, 3 metal/PMMA/Si stacks (**held out**) | 19.6% | 22.3% | 79.4% | 10.2% | 39.2% |
| Chang 2019 TiC enhancement ratio, measured 1.48 (**held out**) | 1.36 | 1.97 | 1.10 | 1.80 | 1.08 |

Sandell stacks (MW m⁻² K⁻¹, measured / model): Au/PMMA/Si 59 / 57, Au/Ti/PMMA/Si 115 / 96, Au/Ni/PMMA/Si 139 / 84

### B. Coherent crystalline interfaces: fast tier vs MACE-AGF (300 K, harmonic)

The AGF tier (`agf.py`, `tbc.py --agf`) is the reference for perfect coherent interfaces. Independent checks: AlN/GaN 281 vs ~300 MW m⁻² K⁻¹ from DFT-NEGF (Polanco & Lindsay, PRB 97, 014303, 2018); AlN/4H-SiC 487 vs ~470 measured (TDTR).

| interface | strain in AGF | **AGF** | fast tier (elastic) | full-dispersion DMM | DMM (Debye) | AMM (Debye) |
|---|---|---|---|---|---|---|
| AlN(0001)/GaN(0001) | -2.9% | **281** | 243 (-14%) | 445 (+59%) | 1376 (+390%) | 3767 (+1241%) |
| AlN(0001)/SiC4H(0001) | +1.1% | **487** | 531 (+9%) | 755 (+55%) | 1849 (+280%) | 2940 (+504%) |
| GaAs(100)/AlAs(100) | +0.7% | **170** | 133 (-21%) | 171 (+1%) | 594 (+250%) | 898 (+429%) |
| GaN(0001)/SiC4H(0001) | +4.1% | **147** | 184 (+25%) | 374 (+154%) | 1291 (+776%) | 1058 (+618%) |
| Si(100)/Ge(100) | -5.2% | **126** | 135 (+8%) | 204 (+62%) | 547 (+335%) | 2086 (+1560%) |
| **mean abs. deviation from AGF** | | | 15% | 66% | 406% | 870% |

### C. Measured conductance of epitaxial semiconductor and nitride/oxide interfaces (300 K, held out)

Multi-frequency TDTR on GaN-on-SiC, GaN-on-AlN and AlN-on-SiC wafers (par.nsf.gov 10498482). * Epitaxial TiN on MgO(001), MgO(111) and sapphire: Costescu, Wall & Cahill, PRB 67, 054302 (2003) found the same conductance on all three; ~700 MW m⁻² K⁻¹ at room temperature as cited by Wilson et al., PRB 91, 115414 (2015) (approximate; TiN electron-phonon term omitted: no g data). Real interfaces contain defects and anharmonic transport; harmonic AGF is a perfect-interface reference.

| interface | measured | fast tier (total) | fast tier (elastic) | AGF (harmonic) | DMM (Debye) | AMM (Debye) |
|---|---|---|---|---|---|---|
| AlN(0001)/GaN(0001) | 490 | 378 | 247 | 281 | 1376 | 3767 |
| AlN(0001)/4H-SiC(0001) | 470 | 657 | 533 | 487 | 1849 | 2940 |
| GaN(0001)/4H-SiC(0001) | 305 | 312 | 189 | 147 | 1291 | 1058 |
| TiN(001)/MgO(001) * | 700 | 511 | 366 | – | 2167 | 3932 |
| TiN(111)/MgO(111) * | 700 | 495 | 337 | – | 2167 | 3932 |
| TiN(111)/sapphire(0001) * | 700 | 482 | 372 | – | 2088 | 3497 |
| **mean abs. error** | | 25% | | 33% | 236% | 461% |

![parity](results/parity.png)

*Parity plot of every validated case: experiments (Li, Cheng, Sandell; all temperatures), measured epitaxial interfaces (section C) and MACE-AGF references (section B), for this model and for the Debye AMM/DMM.*
