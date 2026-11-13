# DiscoveredMaterials — Grand Challenge 05 submission

**Predicting heat transport across solid interfaces**: a two-tier framework
from atomic structure to thermal boundary conductance G(T) and resistance
R_K(T) = 1/G(T), with no per-specimen fitting.

* **Fast tier:** the full-dispersion Landauer model on MACE-OMAT-0 phonons,
  taking seconds per interface.
* **Accurate tier:** an atomistic Green's function with MACE force constants of
  the actual interface, ~1 CPU-hour.

Everything is in [`tbc/`](tbc/):

* [`tbc/README.md`](tbc/README.md): method, all results with figures, and an
  honest list of limitations.
* [`tbc/SUBMISSION.md`](tbc/SUBMISSION.md): the one-page submission summary.
* [`tbc/EXPERIMENTS.md`](tbc/EXPERIMENTS.md): the log of every tested hypothesis.

```bash
cd tbc && pip install -r requirements.txt
python tbc.py AlN/0001 GaN/0001 --T 100:800:50     # any interface
python make_report.py                              # regenerate every quoted number
python -m pytest -q tests/                         # physics invariants
```
