"""Calibrate the two global parameters (alpha, b) of the inelastic channel on
ONE dataset only: Al(111)/Si(111) sample 1 of Li et al. (2022). Everything
else (Al/GaN, Cu/TiC/diamond, Au/PMMA stacks) is a prediction.
Writes params.json."""
import json, sys, os
import numpy as np
from scipy.optimize import minimize
from predict import interface, G_ep, PARAM_FILE, bonding_for
import model as _model
# H15: also fit T_L of the measured-temperature convention f(T) = T/(T+T_L)
FIT_TL = os.environ.get("TBC_FIT_TL", "0") == "1"
LE_N = float(os.environ["TBC_LE_N"]) if os.environ.get("TBC_LE_N") else None
_model.LOCAL_EQ_EXPERIMENT = FIT_TL
_model.LOCAL_EQ_N = LE_N

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


FREEZE_AB = os.environ.get("TBC_FREEZE_AB", "0") == "1"   # H15c: keep alpha, b; fit only T_L
_P0 = json.load(open(PARAM_FILE)) if os.path.exists(PARAM_FILE) else {}


def loss(x):
    if FREEZE_AB:
        P = dict(alpha=_P0["alpha"], b=_P0["b"]); _model.LOCAL_EQ_TL = float(np.exp(x[0]))
    else:
        P = dict(alpha=np.exp(x[0]), b=np.exp(x[1]))
        if FIT_TL:
            _model.LOCAL_EQ_TL = float(np.exp(x[2]))
    c = curve(P, d[:, 0])
    return np.mean(np.log(c / d[:, 1]) ** 2) + W_SLOPE * (slope(d[hi, 0], c[hi]) - slope(d[hi, 0], d[hi, 1])) ** 2


best = None
if FREEZE_AB:
    for t0 in [np.log(t) for t in (500, 650, 800, 1000, 1300)]:
        r = minimize(lambda x: loss(x), [t0], method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-10, maxiter=400))
        if best is None or r.fun < best.fun:
            best = r
    best.x = np.array([np.log(_P0["alpha"]), np.log(_P0["b"]), best.x[0]])
starts_TL = [np.log(t) for t in (150, 300, 600, 1200)] if FIT_TL else [None]
for a0 in ([] if FREEZE_AB else [-1, 0, 1, 2]):
    for b0 in [-1, 0, 1, 2, 3]:
        for t0 in starts_TL:
            x0 = [a0, b0] + ([t0] if t0 is not None else [])
            r = minimize(loss, x0, method="Nelder-Mead",
                         options=dict(xatol=1e-5, fatol=1e-10, maxiter=4000))
            if best is None or r.fun < best.fun:
                best = r
old = json.load(open(PARAM_FILE)) if os.path.exists(PARAM_FILE) else {}
P = dict(alpha=round(float(np.exp(best.x[0])), 4), b=round(float(np.exp(best.x[1])), 4),
         calibrated_on="Li et al. 2022 Al/Si sample 1 (12 points, 82-705 K); objective: log error + %g x (slope >400 K error)^2" % W_SLOPE,
         rms_log_error=round(float(np.sqrt(best.fun)), 4))
if FIT_TL:
    P["local_eq_TL"] = round(float(np.exp(best.x[2])), 1)
    P["local_eq_n"] = LE_N
    P["local_eq_note"] = ("measured G in a partial local-equilibrium convention: fraction f(T) = 1/(1+(T_L/T)^n) of the "
                          "Zeng-Chen correction (H15c); n fixed, T_L fitted on the Al/Si curve" + (" with alpha, b frozen" if FREEZE_AB else ""))
    _model.LOCAL_EQ_TL = P["local_eq_TL"]
    c = curve(P, d[:, 0])
    P["fit_check"] = dict(G300=round(float(np.interp(300, d[:, 0], c)), 1),
                          slope_above_400K=round(float(slope(d[hi, 0], c[hi])), 3), measured_slope=round(float(slope(d[hi, 0], d[hi, 1])), 3))
for k, v in old.items():          # keep every other global parameter (eta, w, ...)
    if not (k in ("local_eq_TL", "local_eq_note", "fit_check") and not FIT_TL):
        P.setdefault(k, v)
json.dump(P, open(PARAM_FILE, "w"), indent=1)
print(json.dumps(P))
