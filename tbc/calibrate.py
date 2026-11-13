"""Calibrate the two global parameters (alpha, b) of the inelastic channel on
ONE dataset only: Al(111)/Si(111) sample 1 of Li et al. (2022). Everything
else (Al/GaN, Cu/TiC/diamond, Au/PMMA stacks) is a prediction.
Writes params.json."""
import json, sys, os
import numpy as np
from scipy.optimize import minimize
from predict import interface, G_ep, PARAM_FILE, bonding_for

ref = json.load(open("data/reference/li2022_fig1.json"))["data"]
d = np.array(ref["Al/Si sample1"])
I = interface("Al", "Si", "111", "111", bonding=bonding_for("Al", "Si", "111", "111"))


def curve(P, Ts):
    return np.array([1 / (1 / G_ep("Al", T) + 1 / I.G_pp(T, P["alpha"], P["b"])) / 1e6 for T in Ts])


# The challenge asks to reproduce "the room-temperature values and the renewed,
# approximately linear increase above about 400 K": the objective is the mean
# squared log error plus a penalty on the slope above 400 K (MW m^-2 K^-2).
W_SLOPE = float(os.environ.get("TBC_SLOPE_WEIGHT", "0"))
hi = d[:, 0] >= 390


def slope(T, g):
    return np.polyfit(T, g, 1)[0]


def loss(x):
    P = dict(alpha=np.exp(x[0]), b=np.exp(x[1]))
    c = curve(P, d[:, 0])
    return np.mean(np.log(c / d[:, 1]) ** 2) + W_SLOPE * (slope(d[hi, 0], c[hi]) - slope(d[hi, 0], d[hi, 1])) ** 2


best = None
for a0 in [-1, 0, 1, 2, 3]:
    for b0 in [-1, 0, 1, 2]:
        r = minimize(loss, [a0, b0], method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-10, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
old = json.load(open(PARAM_FILE)) if os.path.exists(PARAM_FILE) else {}
P = dict(alpha=round(float(np.exp(best.x[0])), 4), b=round(float(np.exp(best.x[1])), 4),
         calibrated_on="Li et al. 2022 Al/Si sample 1 (12 points, 82-705 K); objective: log error + %g x (slope >400 K error)^2" % W_SLOPE,
         rms_log_error=round(float(np.sqrt(best.fun)), 4))
for k in ("spec_efficiency", "spec_efficiency_note"):
    if k in old:
        P[k] = old[k]
json.dump(P, open(PARAM_FILE, "w"), indent=1)
print(json.dumps(P))
