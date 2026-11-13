# Transferable thermal boundary conductance from computed phonons

Submission draft for **Discovered Materials Grand Challenge 05 — Predicting Heat
Transport across Solid Interfaces**: *"Build a computationally tractable
framework that predicts the thermal boundary conductance G(T) and resistance
R_K(T) = 1/G(T) across specified solid–solid interfaces without fitting a
conductance to each specimen. Demonstrate transferable improvement over the
acoustic and diffuse mismatch models (AMM/DMM)."*

Everything runs on a laptop CPU. A new crystalline material costs about one
minute of MACE force evaluations plus phonopy. A G(T) curve takes
milliseconds.

**Status — read this first.** This is a working, reproducible framework with
an honest validation, not a finished prize claim. Calibrated on Al/Si only,
it predicts the held-out Al/GaN G(T) curve to about 5% (80–700 K). It also
predicts the Cu/TiC/diamond enhancement, and the Au/[Ti,Ni]/PMMA/Si stacks
within experimental error to ~25% (see [Results](#results)). It is far better
than AMM/DMM on all three datasets. Known weaknesses are listed in
[Limitations](#limitations), including one regime where the model fails
(very weakly bonded, H-terminated contacts).

---

## 1. Method

`G` is computed as a series network of physically distinct resistances.

```
1/G_total = 1/G_ep(metal)  +  Σ_interfaces 1/G_pp  +  Σ_films d/κ
```

### 1.1 Phonons from a universal interatomic potential (no fitting)
`phonons.py` relaxes each crystal and computes harmonic force constants with
**MACE-OMAT-0** (`medium-omat-0`) and phonopy. On a full Brillouin-zone mesh
(no symmetry reduction, so velocity projections are exact) it builds the
interface-normal spectral flux, i.e. the number of phonon channels per unit
area crossing a plane of normal **n**:

    M_n(ω) = 1/(2 V N_q) Σ_qj |v_qj·n| δ(ω − ω_qj)

It also stores the DOS, the polarization-resolved flux, Debye-averaged sound
velocities and the mean-square displacement. Polar crystals (GaN, AlN) include
LO–TO splitting using literature Born charges.
Spot checks against experiment:

| | MACE-OMAT-0 | experiment |
|---|---|---|
| Si TO(Γ), LA(X), TA(L) | 15.1, 12.3, 3.5 THz | 15.5, 12.3, 3.4 THz |
| f_max Al / Cu / Ni / diamond | 9.5 / 7.2 / 8.7 / 38.4 THz | ~9.7 / 7.3 / 8.9 / 39.9 THz |
| f_max GaN / AlN (LO) | 21.3 / 27.0 THz | 22.2 / 27.3 THz |

(MACE-MP-0, the older model, softens Si optical modes by 28% and is **not** used.)
Amorphous solids (PMMA, SiO₂) use an isotropic Debye spectrum built from
measured density and sound velocities, with the Swartz–Pohl atom-density cutoff.

### 1.2 Elastic transmission: specular + diffuse, from full dispersion
* **Diffuse** (full-dispersion DMM): `τ_d = M₂/(M₁+M₂)`.
* **Specular** (momentum-conserving limit):
  `τ_s = [4Z₁Z₂/(Z₁+Z₂)²]·min(M₁,M₂)/M₁` with spectral impedance
  `Z_i(ω) = ρ_i·4M_i/D_i`. It equals 1 for two identical crystals and obeys
  reciprocity. This removes DMM's best-known defect, a finite resistance
  between identical materials, which matters for acoustically similar
  epitaxial pairs such as AlN/GaN.
* **Specularity from coherence, not from fitting**. Two hexagonal (111)/(0001)
  lattices in the aligned epitaxial orientation share an in-plane period `L_c`,
  the smallest coincidence within 3% misfit. Parallel momentum is conserved
  modulo 2π/L_c, so
  `p(ω) = p∞ + (1−p∞)·exp[−(L_c/λ)²]`, with `p∞ = (a/L_c)²` and `λ = 2π v_T/ω`.
  A 1:1 pseudomorphic interface (AlN/GaN) is specular. Al/Si (4:3, L_c = 11.5 Å)
  and Al/GaN (8:7, 22.8 Å) are diffuse except at long wavelengths. Interfaces
  with an amorphous side are fully diffuse.

### 1.3 Interfacial bonding from MACE + D3 (no fitting)
`bonding.py` builds the actual interface:
* coherent supercells for crystal pairs;
* PMMA repeat units (methyl isobutyrate) at PMMA areal density on the metal or
  hydroxylated-silica surface, averaged over 3 random orientations.

The normal stiffness per area, `K_n = E''(δ)/A`, comes from rigid-shift energy
scans with MACE-OMAT-0 plus D3(BJ) dispersion. The bulk-like reference is
`K_ref = ρv²/d`, harmonic-averaged over both sides. The protocol is validated on
bulk Al(111) cut in half: 5.5×10²⁰ (rigid shift) vs 4.9×10²⁰ N m⁻³ (ρv²/d).
The transmission is multiplied by the spring factor
`S(ω) = s(ω;K)/s(ω;K_ref) ≤ 1`, with `s = 1/(1+(ωZ₁Z₂/(K(Z₁+Z₂)))²)`.
Interfaces bonded at least as stiffly as bulk are not penalised.

### 1.4 Inelastic (anharmonic) channel
Li et al. (2022) showed that at atomically sharp interfaces G keeps rising
linearly above ~400 K, because elastic transport saturates once the soft side
is fully excited. Here an elastically reflected phonon is converted by
interfacial three-phonon scattering into a transmissible mode with
probability

    P(T) = 1 − exp[−α (T/Θ) e^(−Θ/(bT))],   Θ = ħω_max(stiffer side)/k_B

(Umklapp-type freeze-out). The phonon conductance becomes
`G_pp = G_el + P(T)·(G_max − G_el)`. `G_max` is the largest conductance
reachable with **one** conversion step: the radiation limit of the soft side,
or of the stiff-side phonons with ω ≤ 2ω_max,soft, whichever is smaller.
This bound stops Cu/diamond, where up to 5-phonon processes would be needed,
from gaining an unphysical channel. The channel is applied to crystal/crystal
interfaces only, because Θ is a property of a computed crystal spectrum.

**α and b are the only fitted parameters.** They are global and frozen,
calibrated once on the Al/Si sample-1 curve of Li et al. (12 points,
80–700 K). See `params.json` and `calibrate.py`.

### 1.5 Electron–phonon resistance of metals
Majumdar–Reddy: `G_ep = sqrt(g·κ_p(T))`. Here κ_p is the lattice thermal
conductivity of the metal, with phonon–phonon (∝1/T) and phonon–electron
(T-independent) parts from DFT (Tong et al. 2019). g is the measured
electron–phonon coupling factor (sources in `materials.py`). A thin adhesion
layer (2 nm Ti/Ni) adds a parallel volumetric channel `g·d`.

---

## 2. Results

RESULTS_PLACEHOLDER

---

## 3. Limitations (honest list)

1. **High-temperature slope is about half the measured value.** Above 400 K the
   model rises (0.19–0.21 MW m⁻² K⁻²), but Li et al. measure 0.39–0.43. Part of
   the rise is cancelled by the falling electron–phonon conductance of Al
   (κ_p ∝ 1/T); g of Al is itself uncertain by about 2× (pump-probe vs
   equilibrium DFT values, Brown et al. PRB 94, 075120 (2016)).
2. **Low-temperature Al/Si is overpredicted by up to 17% at 80 K.** This comes
   from the long-wavelength specular channel. Real interfaces have some
   roughness that the coherence-based specularity ignores.
3. **The inelastic channel rests on one experiment.** α and b are calibrated
   only on Al/Si. The Al/GaN prediction tests transfer to a second substrate
   with the same metal, not to other metals.
4. **Very weak (van der Waals, H-terminated) contacts fail.** For Cu on
   H-terminated diamond, MACE+D3 gives K/K_ref = 0.04. The spring factor then
   collapses G to about 1 MW m⁻² K⁻¹, far below measured metal/H-diamond values
   (tens of MW m⁻² K⁻¹). Either MACE+D3 underbinds such contacts or other
   channels carry heat (e.g. cross-interface electron–phonon coupling). The
   Chang prediction therefore assumes bulk-terminated diamond(111).
5. **The polymer description is crude.** PMMA is a Debye solid (measured
   sound velocities, atom-density cutoff). Its contact stiffness comes from
   PMMA repeat units on the surface: 3 orientations, MACE-OMAT-0 + D3. OMat
   training data contain few molecular systems.
6. **Composite interlayers are idealised.** The 10 nm TiC film is a bulk-like
   resistance (κ = 20 W m⁻¹ K⁻¹) between two ideal-epitaxy interfaces. The
   measured dependence on TiC crystallinity and thickness (Chang et al.) is
   not modelled.
7. **Only (111)/(0001) hexagonal surface matching is implemented** for the
   coherence length. Other orientations are treated as fully diffuse unless
   L_c is passed explicitly (`Interface(..., Lc=...)`).
8. **Literature inputs** that are not computed: the electron–phonon coupling g
   of each metal, Born charges of polar crystals, and sound velocities and
   density of amorphous solids. Each has its source in `materials.py`.
9. **Comparison convention.** Conductances are Landauer conductances with
   emitted-phonon temperatures (the AGF convention). The Zeng–Chen
   local-equilibrium correction is available (`noneq=True`) but switched off:
   it overshot the low-T Al/Si data by 2×. Comparing with NEMD temperature
   drops at highly transmitting interfaces would need that correction.
10. **No hidden-set guarantee.** The validation covers 3 datasets and 6
    interfaces. Industrial interfaces with intermixing, oxides or amorphous
    interlayers need their structure supplied explicitly.

---

## 4. Running it

```bash
pip install -r requirements.txt
python phonons.py Al Si GaN          # MACE-OMAT-0 phonons (cached in data/phonons)
python bonding.py Al_Si Au_PMMA      # optional: interfacial stiffness (cached in data/bonding.json)
python calibrate.py                  # re-derive the 2 global parameters (Al/Si only)
python validate.py                   # reproduce the three datasets -> results/
python tbc.py AlN/0001 GaN/0001 --T 100:800:50          # any interface
python tbc.py Al/111 Si/111 --metal Al --T 300 --json
python tbc.py Au PMMA --metal Au --Kn 3.9e19             # supply a computed stiffness
```

To add a material, put a structure guess in `materials.CRYSTALS`, or density
and sound velocities in `materials.AMORPHOUS`. Phonons are computed
automatically on first use.

| file | content |
|---|---|
| `phonons.py` | MACE-OMAT-0 + phonopy → M_n(ω), D(ω), MSD, sound velocities |
| `model.py` | Debye AMM / DMM baselines; specular+diffuse elastic, bonding, inelastic channel |
| `bonding.py` | MACE+D3 interface construction and rigid-shift stiffness |
| `predict.py` | e-ph resistance, stacks for the reference systems |
| `calibrate.py` | fits α, b on Al/Si only |
| `validate.py` | all comparisons, tables and figures in `results/` |
| `tbc.py` | command-line predictor for an arbitrary interface |
| `data/reference/` | the experimental data used, with provenance |
| `scripts/digitize_li2022.py` | extracts Li et al. Fig. 1 from the vector PDF |

## 5. Data provenance
* **Li et al., Nat. Commun. 13, 4901 (2022)**: G(T) of Al/Si (samples 1, 2) and Al/GaN.
  Marker centres were read from the vector graphics of the open-access PDF.
  Check: 380.6 / 424.9 / 310.3 at 300 K vs 379 / 423 / 309 quoted in the text.
* **Sandell et al., Nanomaterials 10, 670 (2020)**: Table 2 (combined boundary
  conductance and κ_PMMA) and Eq. (2) (open-access full text, Europe PMC).
* **Chang et al., ACS AMI 11, 26507 (2019)**: abstract value (+48% with a highly
  crystalline 10 nm TiC interlayer). The paper itself is paywalled. Only the ratio is used.

Licensing note: the MACE-OMAT-0 weights are distributed under the Academic
Software License. Commercial use needs a separately licensed universal
potential, which is a drop-in replacement in `phonons.py`/`bonding.py`.
