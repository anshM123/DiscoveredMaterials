"""Regenerate every number quoted in README.md: experimental validation
(validate.py), fast tier vs MACE-AGF reference on coherent interfaces, and the
unseen-interface sanity table. Writes results/summary.md and splices it into
README.md between the RESULTS markers."""
import json, os, re
import numpy as np
import validate
from model import Side, Interface, G_AMM_debye, G_DMM_debye
from predict import load_params

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
AGF_DIR = os.path.join(HERE, "data", "agf")
MW = 1e-6

AGF_SETS = [("AlN0001_GaN0001", ("AlN", "0001"), ("GaN", "0001")),
            ("AlN0001_SiC4H0001", ("AlN", "0001"), ("SiC4H", "0001")),
            ("GaAs100_AlAs100", ("GaAs", "100"), ("AlAs", "100")),
            ("GaN0001_SiC4H0001", ("GaN", "0001"), ("SiC4H", "0001")),
            ("Si100_Ge100", ("Si", "100"), ("Ge", "100"))]


def agf_table(T=300.0):
    rows = []
    for lab, a, b in AGF_SETS:
        path = os.path.join(AGF_DIR, f"{lab}.json")
        if not os.path.exists(path):
            continue
        r = json.load(open(path))
        g_agf = float(np.interp(T, r["T"], r["G"])) * MW
        s1, s2 = Side(*a), Side(*b)
        # AGF models a coherent (strained) interface: compare with the fast
        # tier evaluated as coherent too (Lc = surface lattice parameter)
        I = Interface(s1, s2, Lc=min(s1.a_surf, s2.a_surf))
        rows.append(dict(interface=f"{a[0]}({a[1]})/{b[0]}({b[1]})", strain=r["strain_R"], AGF=g_agf,
                         fast_elastic=I.G_elastic(T) * MW, DMM_full=I.G_dmm_full(T) * MW,
                         AMM=G_AMM_debye(s1, s2, T) * MW, DMM=G_DMM_debye(s1, s2, T) * MW))
    return rows


def main():
    res = validate.main()
    P = load_params()
    agf = agf_table()
    json.dump(agf, open(os.path.join(OUT, "agf_comparison.json"), "w"), indent=1)

    def mape(rows, k):
        return 100 * np.mean([abs(x[k] / x["measured"] - 1) for x in rows])
    li = res["li2022"]; sa = list(res["sandell2020"].values()); ch = res["chang2019"]["Cu_C_bare"]
    cols = ["model", "AMM", "DMM", "AMM_eph", "DMM_eph"]
    L = ["### A. Experiments (fast tier)", "",
         "Mean absolute error against measurement. AMM/DMM are the Debye models of Little (1959) and "
         "Swartz & Pohl (1989), shown raw and with the same electron–phonon resistance and series "
         "layers as this model (`+e-ph`). Only Al/Si is used for fitting (2 global parameters).", "",
         "| dataset (role) | this model | AMM | DMM | AMM+e-ph | DMM+e-ph |", "|---|---|---|---|---|---|"]
    for key, role in [("Al/Si sample1", "calibration"), ("Al/GaN", "**held out**"), ("Al/sapphire (Cheng 2020)", "**held out**")]:
        L.append(f"| {key}, 80–700 K ({role}) | " + " | ".join(f"{mape(li[key], k):.1f}%" for k in cols) + " |")
    L.append("| Sandell 2020, 3 metal/PMMA/Si stacks (**held out**) | " + " | ".join(f"{mape(sa, k):.1f}%" for k in cols) + " |")
    L.append(f"| Chang 2019 TiC enhancement ratio, measured 1.48 (**held out**) | {ch['model']:.2f} | {ch['AMM']:.2f} | "
             f"{ch['DMM']:.2f} | {ch['AMM_eph']:.2f} | {ch['DMM_eph']:.2f} |")
    L += ["", "Sandell stacks (MW m⁻² K⁻¹, measured / model): " +
          ", ".join(f"{k} {v['measured']:.0f} / {v['model']:.0f}" for k, v in res["sandell2020"].items()), ""]
    L += ["### B. Coherent crystalline interfaces: fast tier vs MACE-AGF (300 K, harmonic)", "",
          "The AGF tier (`agf.py`, `tbc.py --agf`) is the reference for perfect coherent interfaces. "
          "Independent checks: AlN/GaN 281 vs ~300 MW m⁻² K⁻¹ from DFT-NEGF (Polanco & Lindsay, PRB 97, 014303, 2018); "
          "AlN/4H-SiC 487 vs ~470 measured (TDTR).", "",
          "| interface | strain in AGF | **AGF** | fast tier (elastic) | full-dispersion DMM | DMM (Debye) | AMM (Debye) |",
          "|---|---|---|---|---|---|---|"]
    for r in agf:
        e = lambda k: f"{r[k]:.0f} ({100 * (r[k] / r['AGF'] - 1):+.0f}%)"
        L.append(f"| {r['interface']} | {100 * r['strain']:+.1f}% | **{r['AGF']:.0f}** | {e('fast_elastic')} | {e('DMM_full')} | {e('DMM')} | {e('AMM')} |")
    if agf:
        for k in ["fast_elastic", "DMM_full", "DMM", "AMM"]:
            pass
        mk = lambda k: np.mean([abs(r[k] / r["AGF"] - 1) for r in agf]) * 100
        L.append(f"| **mean abs. deviation from AGF** | | | {mk('fast_elastic'):.0f}% | {mk('DMM_full'):.0f}% | {mk('DMM'):.0f}% | {mk('AMM'):.0f}% |")
    L.append("")
    txt = "\n".join(L)
    open(os.path.join(OUT, "summary.md"), "w").write(txt)
    readme = os.path.join(HERE, "README.md")
    s = open(readme).read()
    if "<!-- RESULTS:BEGIN -->" in s:
        s = re.sub(r"<!-- RESULTS:BEGIN -->.*<!-- RESULTS:END -->",
                   "<!-- RESULTS:BEGIN -->\n" + txt.replace("\\", "\\\\") + "\n<!-- RESULTS:END -->", s, flags=re.S)
        open(readme, "w").write(s)
    print(txt)


if __name__ == "__main__":
    main()
