# Overview: Grand Challenge 05, heat transport across solid interfaces

**What it is:** a framework that predicts thermal boundary conductance G(T) and
resistance R_K(T) = 1/G(T) for any solid–solid interface, starting from atomic
structure. Nothing is fitted to individual specimens.

**How it works:**
- **Phonons** come from the universal machine-learned potential MACE-OMAT-0.
- **Fast tier**, seconds per curve: a full-dispersion Landauer model with:
  - specular and diffuse transmission;
  - computed interfacial bonding;
  - an anharmonic channel;
  - metal electron terms;
  - thin films and stacks.
- **Accurate tier**, about 1 CPU-hour: an atomistic Green's function on the
  actual relaxed interface.
- **Six global constants** in total, none of them per specimen.

**Headline results:**

| benchmark | this framework | strongest baseline |
|---|---|---|
| Li 2022 Al/Si, 82–705 K (fit set) | 3.6%; 395 vs 381 at 300 K | 18% |
| Li 2022 Al/GaN (not fitted) | 5.7%; 436 vs 425 at 300 K | 13% |
| Cheng 2020 Al/sapphire (not fitted) | 5.8% | 4.6% (Debye AMM + e-ph) |
| Chang 2019 TiC enhancement (1.48) | 1.44 | 1.84 |
| Sandell 2020 metal/PMMA stacks | 13% | 39% |
| 5 pre-registered blind atomistic tests | 5% | 42% |

Units are MW m⁻² K⁻¹. Unless noted, the strongest baseline is full-dispersion
DMM on the same phonons with the same electron–phonon term.

**Reproduce everything:**

```bash
cd tbc && pip install -r requirements.txt
python make_report.py        # regenerates every reported number and figure
python tbc.py AlN/0001 GaN/0001 --T 100:800:50
```

**Where to read more:**
- `tbc/SUBMISSION.md`: one-page summary
- `tbc/README.md`: method, results and limitations
- `tbc/EXPERIMENTS.md`: every hypothesis tested, kept or rejected

## Contributors

- Ansh Mishra
- Aryan Senthilkumar
