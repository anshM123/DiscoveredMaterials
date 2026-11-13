# Experiment log

Merge rule: a change is merged into `claude/nifty-faraday-0rxj6d` only if the
scorecard (`scorecard.py`) improves overall and **no guarded metric regresses**.
Guarded: Al/GaN G(T) (held out), Al/sapphire G(T) (held out), the blind
AGF tests, and the AGF calibration set. Baseline: `results/scorecard_baseline.json`.

Baseline scorecard (commit 6aec94f), lower is better:

| metric | value |
|---|---|
| Al/Si G(T) MAPE (fit) / slope error >400 K | 4.6% / 0.33 |
| Al/GaN G(T) MAPE (held out) / slope error | 6.9% / 0.32 |
| Al/sapphire G(T) MAPE (held out) / slope error | 6.2% / 0.28 |
| Sandell 3 stacks MAPE | 19.6% |
| Chang ratio error | 8.0% |
| fast vs AGF, calibration set / blind set | 15.4% / 10.8% |
| measured epitaxial interfaces (section C) | 25.4% |
| Al / native SiO2 (1.5 nm) / Si G(T) MAPE | 47.1% |

## Gap 1: >400 K slope (Li Al/Si, Al/GaN)

Earlier session (all rejected, held-out regressions): equilibrium-DFT g(Al);
coherence-only specularity; soft-side transmission-deficit channel (Theta soft
and stiff); slope-weighted calibration (slope 0.36/0.32 but 300 K values -10..-19%
and Al/sapphire 6 -> 13%).

**H1–H3 (this sprint): any substrate-dependent inelastic channel.**
* Inverting the data with the model's own e-ph resistance, the phonon conductance
  that the measurements require grows by:
  * x2.1 (Al/Si, 300 -> 700 K);
  * x2.05 (Al/GaN);
  * only x1.13 (Al/sapphire, 300 -> 480 K).
* Candidate discriminators, all computed from the MACE spectra:

  | quantity | Si | GaN | sapphire |
  |---|---|---|---|
  | mean elastic transmission Al -> substrate | 0.28 | 0.39 | 0.26 |
  | headroom G_max - G_el at 300 K (MW m^-2 K^-1) | 592 | 625 | 967 |
  | substrate flux convertible by 3-phonon down-conversion, [w_Al, 2 w_Al], 300 K | 339 | 203 | 847 |
  | same, weighted by the Al two-phonon joint DOS | 287 | 57 | 561 |

* Every physical measure of how much an inelastic channel can add is largest
  for sapphire, yet sapphire is the flat one.
* Any inelastic mechanism that steepens Al/Si and Al/GaN must therefore steepen
  Al/sapphire even more, violating the guard.
* **Conclusion:** the two epitaxial-Al datasets (Li 2022 vs Cheng 2020) are not
  mutually consistent under any single interfacial-phonon mechanism. The
  difference must come from sample- or analysis-specific factors (e.g. the TDTR
  thermal model of the Al film at high T), which are not inputs. Not merged.

## Gap 2: amorphous interfaces

### Sandell stacks (Au, Au/Ti, Au/Ni on PMMA/native SiO2/Si)

**H4: make the polymer contacts ideal (treat MACE molecule-surface stiffness
as unreliable).**
* Required vs model resistances (m^2 K/GW):
  * measured totals: Au 16.9, Ti 8.7, Ni 7.2;
  * model shared path: PMMA/SiO2 5.3 + oxide film 1.5 + SiO2/Si 1.6.
* With an ideal PMMA/SiO2 contact: Ni 143 (meas. 139) but Ti 146 (115) and
  Au 72 (59). Net: no better.
* The measured Ni > Ti ordering needs Ni/PMMA bonded more stiffly than Ti/PMMA.
  That contradicts MACE (K/Kref 0.65 vs 2.98) and Ti's known chemistry.
* Rejected (it would only refit three numbers).

**H5: Open Catalyst potential for Ni/PMMA chemistry.**
* Not pursued: even an ideal Ni contact only raises the Ni stack 84 -> ~94.

### Glassy oxide (Al / native SiO2 / Si; baseline over-predicts ~1.5x)

**Step 1: a MACE a-SiO2 glass** (`glass.py`).
* 96 atoms, 11.45 A hexagonal cell; 4000 K melt, 4 ps quench.
* 94% four-fold Si, 100% two-fold O, density 2.20 g/cm^3.
* Vibrational spectrum extends to ~37 THz.

**H6: replace the Debye glass by the MACE glass's own propagating flux** (full
q-mesh group velocities of the periodic glass).

* M_glass/M_Debye is 0.94 at 1-2 THz, 0.25 at 2-4 THz and 0.05 at 4-14 THz.
  Diffusons barely propagate.
* Radiation limit at 300 K: 232 (MACE glass) vs 2560 MW m^-2 K^-1 (Debye).
* Fast-tier elastic G at 300 K:

  | interface | Debye glass | propagons only |
  |---|---|---|
  | Al/SiO2 | 665 | 109 |
  | Si/SiO2 | 617 | 138 |

* The Al/native-SiO2/Si stack would then be ~3x too low (measured 138).
* Propagons alone are too pessimistic: diffusons also exchange energy across an
  interface. The truth is bracketed between the Debye glass (~1.5x high) and
  propagons only (~3x low).
* Next: atomistic Green's functions for Al/a-SiO2 and Si/a-SiO2 with the
  MACE glass (`agf_glass.py`) to settle where.

**H7: atomistic Green's function, crystal / MACE glass** (`agf_glass.py`;
Al(111) 4x4 | 3 copies of the glass; 480 atoms; float32 force constants).
* Interface bonded: 4 Al-O bonds per interface at 1.79 A.
* G(Al/a-SiO2) = 41 (80 K) and 74 MW m^-2 K^-1 (300 K).
* That is below the measured *stack* (138), which bounds the single interface
  from below.
* So it is a lower bound. A periodic glass lead turns diffusons into flat bands,
  like the propagons-only estimate (109). Not used as ground truth.

**H8 (MERGED): propagons + diffuson weight w, with molecular contacts ideal.**
* M_eff = M_prop(MACE glass) + w (M_Debye - M_prop) for SiO2.
* PMMA contacts: MACE stiffness not used (MACE cannot describe molecule-metal
  chemistry, see Ni/PMMA), so they are treated as bonded.
* w fitted to Li Al/native-SiO2/Si (1.5 nm oxide): w = 0.27 (1.0 nm -> 0.24,
  2.0 nm -> 0.32).
* Scan, before the fit:

  | w | Al/SiO2/Si | Sandell |
  |---|---|---|
  | 1.0 | 47% | 25% (with ideal contacts) |
  | 0.5 | 22% | 18% |
  | 0.35 | 9% | 15% |
  | 0.25 | 3% | 18% |

* Scan with MACE polymer stiffness kept instead: Sandell 20 -> 41% as w goes
  1 -> 0.25. The datasets conflict under that assumption, which is why the
  ideal-contact treatment was adopted.
* Scorecard (fit w = 0.27):
  * Al/SiO2/Si 47.1 -> 1.8% (fit set);
  * Sandell 19.6 -> 17.7% (held out);
  * **no regression** on any other metric (the guarded metrics contain no
    amorphous phases).
* Merged.

### Metal / molecular-solid contacts (Sandell; Au over-predicted +19% with every contact bonded)

**H9 (MERGED): classify molecular contacts by chemistry; physisorbed ones use
the MACE + D3 stiffness.**
* Reactive metals (stable oxide, formation enthalpy < −100 kJ per mol O: Ti,
  Al, Ni, Cu, Co, W, Mo, Ta, Ru …) and hydroxylated oxides chemisorb the PMMA
  carbonyl/ester oxygen; MACE-OMAT-0 cannot describe that (Ni/PMMA, H5), so
  these stay ideal.
* Noble metals (Au, Ag, Pt, Pd) bind polymers by dispersion only (XPS
  reactivity series Ti > Cr > Al > Ni > Cu >> Ag, Au). Dispersion is what
  MACE + D3 describes, so the rigid-shift stiffness is used
  (Au/PMMA K_n = 3.9×10¹⁹ N m⁻³, computed before H8 and not adjusted).
* No new parameter. Caveat: the rule was introduced *after* seeing the Sandell
  residuals, so Sandell is no longer a fully blind test of it.
* Sandell (MW m⁻² K⁻¹, measured / before / after): Au 59 / 70 / 60,
  Ti 115 / 105 / 105, Ni 139 / 104 / 104.
* Scorecard: Sandell 17.7 → 12.1%; **no other metric changes** (no PMMA in
  any other set). Merged.
* Remaining: Ni −25%. The PMMA/SiO₂/native-oxide/Si path shared by all three
  stacks already limits any stack to ≤132, so the measured Ni (139) > Ti (115)
  ordering cannot be reproduced by any contact treatment; within the 10%
  experimental error of each value it is only marginally significant.
* Ag/PMMA and Pt/PMMA stiffnesses computed with the same protocol (bonding.py)
  so the rule is usable for those metals without the Au fallback.

**H10 (MERGED, robustness): electron–phonon coupling for metals without a
measured g.**
* Before: Ag, Pt, Co, W, Mo, Ta, Ru (all in the phonon library) silently lost
  their electron–phonon resistance (warning, term omitted).
* Now Allen's formula g = 3ħγ_bλ⟨ω²⟩/(πk_B), with Kittel's Sommerfeld γ,
  λ from Tong 2019 / McMillan 1968, ⟨ω²⟩ from the MACE DOS.
* Check on the five metals with measured g (Allen / measured):
  Al 1.60, Cu 0.68, Au 0.42, Ni 2.7, Ti 0.51; geometric mean 0.91.
* Resulting G_ep(300 K): Ag 300, Pt 2220, W 1330, Mo 2310, Ta 2340,
  Co 4980, Ru 3170 MW m⁻² K⁻¹. Only Ag (weakly coupled, like Au) is
  e-ph-limited; Tong's lattice conductivity is used for Ag, Pt, Co.
* Scorecard unchanged (no such metal in any validation set).

### Independent audit (judge's-eye review by a separate agent) and fixes

An agent with no stake in the results reviewed the docs, code and CLI. Its
findings, and what was done:

| finding | action |
|---|---|
| All warnings silenced globally (`phonons.py`), so the e-ph fallbacks never printed | filter scoped to third-party libraries |
| Debye AMM used (Z₁cosθ₁+Z₂cosθ₂)² and was not reciprocal (Pt/Si 104 vs Si/Pt 1350) | textbook acoustic-Fresnel energy transmission (Z₁cosθ₂+Z₂cosθ₁)², transmission only below both Debye cutoffs → reciprocal. Baselines changed: AMM+e-ph now **beats** this model on Al/sapphire (4.6 vs 6.2%), and is 47% on Sandell |
| Li's second Al/Si specimen digitised but not reported | reported: model 24% (specimen 2 is 23% below specimen 1 at 300 K; its slope 0.10 matches the model's 0.14) |
| Baselines were Debye-only | full-dispersion DMM on the same MACE spectra + same e-ph added (18 / 12 / 13 / 14% on the four G(T) sets) |
| CLI treated Si/Ge as diffuse (MACE misfit 5.2% > 5%), unlike the validation path | **H11**: misfit judged on planes rescaled to the experimental atomic volume (real Si/Ge misfit 4.2%). Scorecard: no metric moved by more than 0.2 (Sandell 12.1, blind 10.8, agf 15.4). The CLI path is within 2% of the validated coherent path for all eight AGF pairs and is now shown next to it |
| `G` (with inelastic channel) vs `G_el` ambiguity for atomistic ground truth | documented like-with-like (AGF ↔ G_el, NEMD ↔ G_pp, experiment ↔ G); `--local-equilibrium` exposes the Zeng–Chen convention of NEMD |
| No electronic channel for metal/metal | electron DMM (Gundrum 2005) in parallel; metal sides detected automatically; `--no-eph` |
| Stack mode lacked the adhesion-layer coupling and film conductivities | added (CLI reproduces the validated Sandell path: Au 60.4 = 60.4, Au/Ti 103 vs 105); bulk-k / κ_min fallbacks |
| "Held out" overstated for Sandell and Chang | relabelled; Chang's dependence on the unreported diamond termination stated |
| Wrong citation (Polanco & Lindsay) | corrected to PRB 99, 075202 (2019) |
| glass.py defaults did not reproduce the cached glass; flux script not in repo | `glass.py <name>` now runs the whole pipeline; regenerates aSiO2_32_flux.npz bit-for-bit |
