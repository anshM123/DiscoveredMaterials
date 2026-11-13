"""Calibrate the glass diffuson coupling weight W_DIFFUSON on ONE dataset:
Li et al. 2022 Al / native SiO2 / Si G(T) (8 points, 80-350 K), with the
native-oxide thickness taken as 1.5 nm (not reported; 1.0 and 2.0 nm give the
sensitivity). Sandell et al. (PMMA on native SiO2/Si) remains held out.
Writes 'w_diffuson' into params.json."""
import json, os, sys
import numpy as np
import model
from model import Side, Interface
from predict import G_ep, film_k, load_params, PARAM_FILE

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.array(json.load(open(os.path.join(HERE, "data/reference/li2022_fig1.json")))["data"]["Al/SiO2(native)/Si"])
P = load_params()


def curve(w, t_ox):
    model.W_DIFFUSON = w
    I1 = Interface(Side("Al", "111"), Side("SiO2")); I2 = Interface(Side("SiO2"), Side("Si", "100"))
    return np.array([1 / (1 / G_ep("Al", T) + 1 / I1.G_pp(T, P["alpha"], P["b"]) + t_ox / film_k("SiO2", T)
                          + 1 / I2.G_pp(T, P["alpha"], P["b"])) / 1e6 for T in d[:, 0]])


ws = np.linspace(0.05, 1.0, 96)
res = {}
for t in (1.0e-9, 1.5e-9, 2.0e-9):
    err = [np.mean(np.log(curve(w, t) / d[:, 1]) ** 2) for w in ws]
    wb = float(ws[int(np.argmin(err))])
    c = curve(wb, t)
    res[t] = wb
    print(f"oxide {t * 1e9:.1f} nm: w_diffuson = {wb:.2f}, MAPE {100 * np.mean(abs(c / d[:, 1] - 1)):.1f}%")
Pj = json.load(open(PARAM_FILE))
Pj["w_diffuson"] = round(res[1.5e-9], 3)
Pj["w_diffuson_note"] = ("fit to Li 2022 Al/native-SiO2/Si G(T) assuming 1.5 nm oxide (1.0 nm -> %.2f, 2.0 nm -> %.2f); "
                         "Sandell stacks held out" % (res[1.0e-9], res[2.0e-9]))
json.dump(Pj, open(PARAM_FILE, "w"), indent=1)
print(json.dumps({k: Pj[k] for k in ("w_diffuson", "w_diffuson_note")}))
