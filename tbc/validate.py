"""Reproduce the three published datasets of the challenge and compare with
the acoustic and diffuse mismatch models. Writes results/validation.json,
results/validation.md and figures in results/.

Baselines: Debye AMM (Little 1959) and Debye DMM (Swartz & Pohl 1989) for the
phonon interface. They are reported raw (as usually quoted) and, for a fair
comparison, with the same metal electron-phonon resistance and the same series
layers as the model ("+e-ph").
"""
import json, os
import numpy as np
from predict import (load_params, interface, G_ep, bonding_for, sandell_stack,
                     chang_ratio, baselines)
from model import Side, G_AMM_debye, G_DMM_debye
from materials import FILM_K

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
ref_li = json.load(open(os.path.join(HERE, "data/reference/li2022_fig1.json")))["data"]
ref_ch = json.load(open(os.path.join(HERE, "data/reference/chang2019.json")))
ref_sa = json.load(open(os.path.join(HERE, "data/reference/sandell2020.json")))["systems"]
P = load_params()
MW = 1e-6


def series(*G):
    return 1.0 / sum(1.0 / g for g in G)


def li_block():
    out = {}
    for key, sub, nsub in [("Al/Si sample1", "Si", "111"), ("Al/GaN", "GaN", "0001")]:
        bond = bonding_for("Al", sub, "111", nsub)
        I = interface("Al", sub, "111", nsub, bonding=bond)
        s1, s2 = Side("Al", "111"), Side(sub, nsub)
        rows = []
        for T, Gm in ref_li[key]:
            gep = G_ep("Al", T)
            amm, dmm = G_AMM_debye(s1, s2, T), G_DMM_debye(s1, s2, T)
            rows.append(dict(T=T, measured=Gm,
                             model=series(gep, I.G_pp(T, P["alpha"], P["b"])) * MW,
                             model_elastic_only=series(gep, I.G_pp(T, 0.0)) * MW,
                             AMM=amm * MW, DMM=dmm * MW,
                             AMM_eph=series(gep, amm) * MW, DMM_eph=series(gep, dmm) * MW))
        out[key] = rows
    return out


def err(rows, k):
    m = np.array([r["measured"] for r in rows]); p = np.array([r[k] for r in rows])
    return dict(MAPE=float(np.mean(np.abs(p / m - 1)) * 100), max_abs_pct=float(np.max(np.abs(p / m - 1)) * 100))


def slope_above(rows, k, Tmin=400):
    T = np.array([r["T"] for r in rows]); g = np.array([r[k] for r in rows])
    s = T >= Tmin
    return float(np.polyfit(T[s], g[s], 1)[0])


def chang_block(cu_c_key):
    bond = {("Cu", "C"): bonding_for("Cu", "C", "111", "111", key=cu_c_key),
            ("Cu", "TiC"): bonding_for("Cu", "TiC", "111", "111"),
            ("TiC", "C"): bonding_for("TiC", "C", "111", "111")}
    r, parts = chang_ratio(P=P, bond=bond)
    res = dict(model=r, parts={k: v * MW for k, v in parts.items()})
    T = 300.0
    gep = G_ep("Cu", T)
    S = {n: Side(n, "111") for n in ("Cu", "C", "TiC")}
    for lab, f in [("AMM", G_AMM_debye), ("DMM", G_DMM_debye)]:
        g0 = f(S["Cu"], S["C"], T)
        g1 = series(f(S["Cu"], S["TiC"], T), f(S["TiC"], S["C"], T), FILM_K["TiC"] / 10e-9)
        res[lab] = g1 / g0
        res[lab + "_eph"] = series(gep, g1) / series(gep, g0)
    return res


def sandell_block():
    bond = {k: bonding_for(*k) for k in [("Au", "PMMA"), ("Ti", "PMMA"), ("Ni", "PMMA"), ("PMMA", "SiO2")]}
    out = {}
    T = 300.0
    for name, spec in ref_sa.items():
        ad = spec["adhesion"]
        G, parts = sandell_stack(ad, T=T, P=P, bond=bond)
        Ge, _ = sandell_stack(ad, T=T, P=dict(alpha=0.0, b=1.0), bond=bond)
        Gnb, _ = sandell_stack(ad, T=T, P=P, bond={})
        import model as _m
        _m.INELASTIC_AMORPHOUS = True
        Gia, _ = sandell_stack(ad, T=T, P=P, bond=bond)
        _m.INELASTIC_AMORPHOUS = False
        top = ad or "Au"
        b = {}
        for lab, f in [("AMM", G_AMM_debye), ("DMM", G_DMM_debye)]:
            R = sum(1 / f(Side(x), Side(y), T) for x, y in [(top, "PMMA"), ("PMMA", "SiO2"), ("SiO2", "Si")])
            R += 2.05e-9 / FILM_K["SiO2"]
            Rf = R + (spec["adhesion_nm"] * 1e-9 / FILM_K[ad] if ad else 0.0)
            gep = G_ep("Au", T) + (G_ep(ad, T, d=spec["adhesion_nm"] * 1e-9) if ad else 0.0)
            b[lab] = 1 / R * MW
            b[lab + "_eph"] = 1 / (Rf + 1 / gep) * MW
        out[name] = dict(measured=spec["G"], model=G * MW, model_elastic_only=Ge * MW,
                         model_no_bonding=Gnb * MW, model_inelastic_on_amorphous=Gia * MW, parts={k: v * MW for k, v in parts.items()}, **b)
        # thickness dependence G_eff(t) = 1/(t/k + 1/G), k = 0.20 W/mK
        t = np.linspace(1, 15, 15) * 1e-9
        out[name]["Geff_t"] = dict(t_nm=list(t * 1e9),
                                   measured_fit=list(1 / (t / spec["kPMMA"] + 1 / (spec["G"] * 1e6)) * MW),
                                   model=list(1 / (t / FILM_K["PMMA"] + 1 / G) * MW))
    return out


def main():
    res = dict(params=P, li2022=li_block(), sandell2020=sandell_block(),
               chang2019={k: chang_block(k) for k in ("Cu_C_H", "Cu_C_bare")})
    json.dump(res, open(os.path.join(OUT, "validation.json"), "w"), indent=1)
    write_md(res)
    plots(res)
    return res


def write_md(res):
    L = ["# Validation against the three published datasets", "",
         f"Global parameters (calibrated on Al/Si sample 1 only): alpha = {P['alpha']}, b = {P['b']}", ""]
    for key, rows in res["li2022"].items():
        tag = "calibration set" if "Si" in key else "PREDICTION (not fitted)"
        L += [f"## Li et al. 2022 - {key} ({tag})", "",
              "| T (K) | measured | model | elastic only | AMM | DMM | AMM+e-ph | DMM+e-ph |", "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            L.append("| {T:.0f} | {measured:.0f} | **{model:.0f}** | {model_elastic_only:.0f} | {AMM:.0f} | {DMM:.0f} | {AMM_eph:.0f} | {DMM_eph:.0f} |".format(**r))
        L.append("")
        L.append("| metric | model | elastic only | AMM | DMM | AMM+e-ph | DMM+e-ph |")
        L.append("|---|---|---|---|---|---|---|")
        L.append("| mean abs. error (%) | " + " | ".join(f"{err(rows, k)['MAPE']:.1f}" for k in ["model", "model_elastic_only", "AMM", "DMM", "AMM_eph", "DMM_eph"]) + " |")
        L.append("| max abs. error (%) | " + " | ".join(f"{err(rows, k)['max_abs_pct']:.1f}" for k in ["model", "model_elastic_only", "AMM", "DMM", "AMM_eph", "DMM_eph"]) + " |")
        L.append("| slope above 400 K (MW m^-2 K^-2), measured " + f"{slope_above(rows, 'measured'):.2f} | " + " | ".join(f"{slope_above(rows, k):.2f}" for k in ["model", "model_elastic_only", "AMM", "DMM", "AMM_eph", "DMM_eph"]) + " |")
        L.append("")
    L += ["## Sandell et al. 2020 - combined metal/PMMA + PMMA/Si boundary conductance at 300 K (PREDICTION)", "",
          "| stack | measured | model | variant: w/o bonding factor | variant: inelastic channel also on amorphous | AMM | DMM | AMM+e-ph | DMM+e-ph |", "|---|---|---|---|---|---|---|---|---|"]
    for k, r in res["sandell2020"].items():
        L.append(f"| {k} | {r['measured']:.0f} | **{r['model']:.0f}** | {r['model_no_bonding']:.0f} | {r['model_inelastic_on_amorphous']:.0f} | {r['AMM']:.0f} | {r['DMM']:.0f} | {r['AMM_eph']:.0f} | {r['DMM_eph']:.0f} |")
    L.append("")
    L += ["## Chang et al. 2019 - G(Cu/TiC 10 nm/diamond) / G(Cu/diamond) at 300 K (PREDICTION), measured 1.48", "",
          "| diamond surface in model | model | AMM | DMM | AMM+e-ph | DMM+e-ph |", "|---|---|---|---|---|---|"]
    for k, r in res["chang2019"].items():
        L.append(f"| {k} | **{r['model']:.2f}** | {r['AMM']:.2f} | {r['DMM']:.2f} | {r['AMM_eph']:.2f} | {r['DMM_eph']:.2f} |")
    L.append("")
    open(os.path.join(OUT, "validation.md"), "w").write("\n".join(L))


def plots(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, (key, rows) in zip(axs, res["li2022"].items()):
        T = [r["T"] for r in rows]
        ax.plot(T, [r["measured"] for r in rows], "o", color="k", label="measured (Li 2022)")
        Tf = np.linspace(80, 720, 65)
        sub, nsub = ("Si", "111") if "Si" in key else ("GaN", "0001")
        I = interface("Al", sub, "111", nsub, bonding=bonding_for("Al", sub, "111", nsub))
        s1, s2 = Side("Al", "111"), Side(sub, nsub)
        ax.plot(Tf, [series(G_ep("Al", t), I.G_pp(t, P["alpha"], P["b"])) * MW for t in Tf], "-", color="C3", lw=2, label="this model")
        ax.plot(Tf, [series(G_ep("Al", t), I.G_pp(t, 0)) * MW for t in Tf], "--", color="C3", lw=1, label="model, elastic only")
        ax.plot(Tf, [G_DMM_debye(s1, s2, t) * MW for t in Tf], ":", color="C0", label="DMM (Debye)")
        ax.plot(Tf, [G_AMM_debye(s1, s2, t) * MW for t in Tf], "-.", color="C2", label="AMM (Debye)")
        ax.plot(Tf, [series(G_ep("Al", t), G_DMM_debye(s1, s2, t)) * MW for t in Tf], ":", color="C0", alpha=0.5, label="DMM + e-ph")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("T (K)"); ax.set_ylabel("G (MW m$^{-2}$ K$^{-1}$)")
        ax.set_title(f"{key}" + (" (calibration)" if "Si" in key else " (prediction)"))
        ax.set_ylim(80, 2500)
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "li2022.png"), dpi=150)
    fig, ax = plt.subplots(figsize=(5, 4))
    for i, (k, r) in enumerate(res["sandell2020"].items()):
        g = r["Geff_t"]
        ax.plot(g["t_nm"], g["measured_fit"], "-", color=f"C{i}", alpha=0.5, lw=4, label=f"{k} measured (fit)")
        ax.plot(g["t_nm"], g["model"], "--", color=f"C{i}", label=f"{k} model")
    ax.set_xlabel("PMMA thickness (nm)"); ax.set_ylabel("G$_{eff}$ (MW m$^{-2}$ K$^{-1}$)")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(OUT, "sandell2020.png"), dpi=150)


if __name__ == "__main__":
    main()
    print(open(os.path.join(OUT, "validation.md")).read())
