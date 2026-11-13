"""Calibrate the specular efficiencies of the fast tier against the MACE-AGF
references for coherent crystalline interfaces, with leave-one-out errors.

Two global constants (model.spec_efficiency): eta for all pairs, eta_hetero for
heterovalent III-V/group-IV pairs. Calibration set: every valid AGF reference
computed so far that is not in a pending blind round, i.e. the original five
(make_report.AGF_SETS) plus the nine first-round blind pairs (AGF_BLIND, whose
results were already known when this two-class form was chosen). The fresh
blind rounds (make_report.AGF_BLIND2) are never used here.
Writes 'spec_efficiency' and 'spec_efficiency_hetero' into params.json (the
inelastic parameters must be re-calibrated afterwards with calibrate.py)."""
import json, os
import numpy as np
import model
from model import Side, Interface
from materials import heterovalent
from make_report import AGF_SETS, AGF_BLIND, AGF_DIR
from predict import PARAM_FILE

T = 300.0
data = []
for lab, a, b in AGF_SETS + AGF_BLIND:
    p = os.path.join(AGF_DIR, f"{lab}.json")
    if os.path.exists(p):
        r = json.load(open(p))
        if r.get("valid", True):
            data.append((lab, a, b, float(np.interp(T, r["T"], r["G"])), "III-V" in heterovalent(a[0], b[0])))

etas = np.linspace(0.2, 1.2, 101)


def fast(a, b, eta):
    model.SPEC_EFFICIENCY = eta; model.SPEC_EFFICIENCY_HETERO = None
    s1, s2 = Side(*a), Side(*b)
    return Interface(s1, s2, Lc=min(s1.a_surf, s2.a_surf)).G_elastic(T)


pred = np.array([[fast(a, b, e) for e in etas] for _, a, b, _, _ in data])
ref = np.array([d[3] for d in data])
err = np.log(pred / ref[:, None]) ** 2
out, notes = {}, []
for cls, key in ((False, "spec_efficiency"), (True, "spec_efficiency_hetero")):
    idx = [i for i, d in enumerate(data) if d[4] == cls]
    best = etas[np.argmin(err[idx].mean(axis=0))]
    loo = []
    for i in idx:
        rest = [j for j in idx if j != i]
        e_i = etas[np.argmin(err[rest].mean(axis=0))]
        loo.append(100 * (pred[i, np.argmin(abs(etas - e_i))] / ref[i] - 1))
    dev = [100 * (pred[i, np.argmin(abs(etas - best))] / ref[i] - 1) for i in idx]
    print(f"{key} = {best:.2f} on {len(idx)} pairs; deviations " + ", ".join(f"{data[i][0]} {d:+.0f}%" for i, d in zip(idx, dev)))
    print(f"   leave-one-out: " + ", ".join(f"{data[i][0]} {x:+.0f}%" for i, x in zip(idx, loo)) + f"; mean abs {np.mean(np.abs(loo)):.0f}%")
    out[key] = round(float(best), 3)
    notes.append(f"{key}: {len(idx)} pairs, LOO mean abs dev {np.mean(np.abs(loo)):.0f}%")
P = json.load(open(PARAM_FILE)); P.update(out)
P["spec_efficiency_note"] = ("fit to MACE-AGF at 300 K, two valence classes (III-V/IV vs all other pairs); " + "; ".join(notes))
json.dump(P, open(PARAM_FILE, "w"), indent=1)
print(json.dumps(out))
