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
