"""Calibrate the global specular efficiency (eta) of the fast tier against the
MACE-AGF references for coherent crystalline interfaces, with leave-one-out
errors. Writes 'spec_efficiency' into params.json (the inelastic parameters
must be re-calibrated afterwards with calibrate.py)."""
import json, os
import numpy as np
import model
from model import Side, Interface
from make_report import AGF_SETS, AGF_DIR
from predict import PARAM_FILE

T = 300.0
data = []
for lab, a, b in AGF_SETS:
    p = os.path.join(AGF_DIR, f"{lab}.json")
    if os.path.exists(p):
        r = json.load(open(p))
        data.append((lab, a, b, float(np.interp(T, r["T"], r["G"]))))


def fast(a, b, eta):
    model.SPEC_EFFICIENCY = eta
    s1, s2 = Side(*a), Side(*b)
    return Interface(s1, s2, Lc=min(s1.a_surf, s2.a_surf)).G_elastic(T)


etas = np.linspace(0.3, 1.0, 71)
pred = np.array([[fast(a, b, e) for e in etas] for _, a, b, _ in data])
ref = np.array([d[3] for d in data])
err = np.log(pred / ref[:, None]) ** 2
best = etas[np.argmin(err.mean(axis=0))]
print("data:", [(d[0], round(d[3] / 1e6)) for d in data])
print(f"eta (all) = {best:.2f}; deviations: " + ", ".join(f"{d[0]} {100 * (pred[i, np.argmin(abs(etas - best))] / ref[i] - 1):+.0f}%" for i, d in enumerate(data)))
loo = []
for i in range(len(data)):
    m = np.ones(len(data), bool); m[i] = False
    e_i = etas[np.argmin(err[m].mean(axis=0))]
    loo.append(100 * (pred[i, np.argmin(abs(etas - e_i))] / ref[i] - 1))
    print(f"  leave-out {data[i][0]}: eta={e_i:.2f}  predicted deviation {loo[-1]:+.0f}%")
print(f"LOO mean abs deviation {np.mean(np.abs(loo)):.0f}%  (eta=1: {np.mean(np.abs(pred[:, -1] / ref - 1)) * 100:.0f}%)")
P = json.load(open(PARAM_FILE)); P["spec_efficiency"] = round(float(best), 3)
P["spec_efficiency_note"] = f"fit to MACE-AGF on {len(data)} coherent pairs at 300 K; LOO mean abs dev {np.mean(np.abs(loo)):.0f}%"
json.dump(P, open(PARAM_FILE, "w"), indent=1)
