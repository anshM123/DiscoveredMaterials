"""Scorecard used for the merge rule of experimental changes.

score() returns every validation metric in one dict (lower is better unless
noted). `python scorecard.py` prints it; `compare(a, b)` lists regressions.
Guarded (held-out / blind) metrics must not regress for a change to be merged.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

GUARDED = ["AlGaN_mape", "Alsapph_mape", "blind_mad", "agf_mad"]


def _slope(T, g, Tmin=390):
    T = np.asarray(T); g = np.asarray(g); s = T >= Tmin
    return float(np.polyfit(T[s], g[s], 1)[0])


def score(P=None, verbose=False):
    import validate, make_report
    from predict import load_params, G_ep, film_k
    from model import Side, Interface
    P = P or load_params()
    validate.P.clear(); validate.P.update(P)
    out = {}
    li = validate.li_block()
    for key, tag in [("Al/Si sample1", "AlSi"), ("Al/Si sample2", "AlSi2"), ("Al/GaN", "AlGaN"), ("Al/sapphire (Cheng 2020)", "Alsapph")]:
        r = li[key]
        m = np.array([x["measured"] for x in r]); p = np.array([x["model"] for x in r]); T = [x["T"] for x in r]
        out[f"{tag}_mape"] = float(np.mean(np.abs(p / m - 1)) * 100)
        out[f"{tag}_slope_err"] = abs(_slope(T, p) - _slope(T, m))
        i300 = int(np.argmin(np.abs(np.array(T) - 300)))
        out[f"{tag}_300K_err"] = float(abs(p[i300] / m[i300] - 1) * 100)
    sa = validate.sandell_block()
    out["Sandell_mape"] = float(np.mean([abs(v["model"] / v["measured"] - 1) for v in sa.values()]) * 100)
    ch = validate.chang_block("Cu_C_bare")
    out["Chang_err"] = float(abs(ch["model"] / 1.48 - 1) * 100)
    agf = make_report.agf_table()
    out["agf_mad"] = float(np.mean([abs(r["fast_elastic"] / r["AGF"] - 1) for r in agf]) * 100)
    blind = make_report.agf_table(sets=make_report.AGF_BLIND)
    out["blind_mad"] = float(np.mean([abs(r["fast_elastic"] / r["AGF"] - 1) for r in blind]) * 100)
    b2 = make_report.agf_table(sets=make_report.AGF_BLIND2)
    if b2:
        out["blind2_mad"] = float(np.mean([abs(r["fast_elastic"] / r["AGF"] - 1) for r in b2]) * 100)
    # section C (measured epitaxial), total conductance
    meas = [(("AlN", "0001"), ("GaN", "0001"), 490), (("AlN", "0001"), ("SiC4H", "0001"), 470),
            (("GaN", "0001"), ("SiC4H", "0001"), 305), (("TiN", "100"), ("MgO", "100"), 700),
            (("TiN", "111"), ("MgO", "111"), 700), (("TiN", "111"), ("Al2O3", "0001"), 700)]
    errs = []
    for a, b, gm in meas:
        I = Interface(Side(*a), Side(*b))
        errs.append(abs(I.G_pp(300.0, P["alpha"], P["b"]) / 1e6 / gm - 1))
    out["secC_mape"] = float(np.mean(errs) * 100)
    # section D: Al / native SiO2 (1.5 nm) / Si
    d = np.array(json.load(open(os.path.join(HERE, "data/reference/li2022_fig1.json")))["data"]["Al/SiO2(native)/Si"])
    I1 = Interface(Side("Al", "111"), Side("SiO2")); I2 = Interface(Side("SiO2"), Side("Si", "100"))
    c = np.array([1 / (1 / G_ep("Al", T) + 1 / I1.G_pp(T, P["alpha"], P["b"]) + 1.5e-9 / film_k("SiO2", T)
                       + 1 / I2.G_pp(T, P["alpha"], P["b"])) / 1e6 for T in d[:, 0]])
    out["AlSiO2Si_mape"] = float(np.mean(np.abs(c / d[:, 1] - 1)) * 100)
    if verbose:
        print(json.dumps({k: round(v, 3) for k, v in out.items()}, indent=1))
    return out


def compare(base, new, tol=0.2):
    """Return (improved, regressed) metric lists; tol in metric units."""
    imp, reg = [], []
    for k in base:
        if new[k] < base[k] - tol:
            imp.append((k, round(base[k], 2), round(new[k], 2)))
        elif new[k] > base[k] + tol:
            reg.append((k, round(base[k], 2), round(new[k], 2)))
    return imp, reg


if __name__ == "__main__":
    s = score(verbose=True)
    if len(sys.argv) > 1 and sys.argv[1] == "--save":
        json.dump(s, open(os.path.join(HERE, "results", "scorecard_baseline.json"), "w"), indent=1)
