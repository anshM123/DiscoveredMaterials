"""Calibrate the two global parameters (alpha, b) of the inelastic channel on
ONE dataset only: Al(111)/Si(111) sample 1 of Li et al. (2022). Everything
else (Al/GaN, Cu/TiC/diamond, Au/PMMA stacks) is a prediction.
Writes params.json."""
import json, sys
import numpy as np
from scipy.optimize import minimize
from predict import interface, G_ep, PARAM_FILE, bonding_for

ref = json.load(open("data/reference/li2022_fig1.json"))["data"]
d = np.array(ref["Al/Si sample1"])
I = interface("Al", "Si", "111", "111", bonding=bonding_for("Al", "Si", "111", "111"))


def curve(P, Ts):
    return np.array([1 / (1 / G_ep("Al", T) + 1 / I.G_pp(T, P["alpha"], P["b"])) / 1e6 for T in Ts])


def loss(x):
    P = dict(alpha=np.exp(x[0]), b=np.exp(x[1]))
    return np.mean(np.log(curve(P, d[:, 0]) / d[:, 1]) ** 2)


best = None
for a0 in [-1, 0, 1, 2, 3]:
    for b0 in [-1, 0, 1]:
        r = minimize(loss, [a0, b0], method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-10, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
P = dict(alpha=round(float(np.exp(best.x[0])), 4), b=round(float(np.exp(best.x[1])), 4),
         calibrated_on="Li et al. 2022 Al/Si sample 1 (12 points, 82-705 K)",
         rms_log_error=round(float(np.sqrt(best.fun)), 4))
json.dump(P, open(PARAM_FILE, "w"), indent=1)
print(json.dumps(P))
