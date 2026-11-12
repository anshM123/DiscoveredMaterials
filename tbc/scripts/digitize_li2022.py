"""Re-create data/reference/li2022_fig1.json from the open-access PDF of
Li et al., Nat. Commun. 13, 4901 (2022) (e.g. https://www.osti.gov/pages/servlets/purl/2470902).

The figure is vector graphics: marker centres are read from the PDF drawing
commands and mapped through log-log axes calibrated on the tick labels.
Check: at 300 K this gives 380.6 / 424.9 / 310.3 MW m^-2 K^-1 vs the values
379 / 423 / 309 quoted in the paper text (Al/Si S1, Al/GaN, Al/Si S2).

usage: python scripts/digitize_li2022.py path/to/li2022.pdf
"""
import json, sys
import numpy as np
import pymupdf

pdf = sys.argv[1]
p = pymupdf.open(pdf)[2]
words = [w for w in p.get_text("words") if w[1] < 245]


def axis(x0, x1):
    xs, ys = {}, {}
    for w in words:
        cx, cy = (w[0] + w[2]) / 2, (w[1] + w[3]) / 2
        if x0 < cx < x1:
            if abs(cy - 236.8) < 1:
                xs[float(w[4])] = cx
            elif w[2] < x0 + 40 and w[4] in ("40", "100", "300", "1000"):
                ys[float(w[4])] = cy
    kx = np.polyfit(np.log10(list(xs)), list(xs.values()), 1)
    ky = np.polyfit(np.log10(list(ys)), list(ys.values()), 1)
    return kx, ky


out = {}
for name, (x0, x1), col in [("Al/Si sample1", (120, 300), (0.81, 0.16, 0.1)),
                            ("Al/Si sample2", (120, 300), (0.27, 0.5, 0.73)),
                            ("Al/GaN", (300, 490), (0.81, 0.16, 0.1))]:
    kx, ky = axis(x0, x1)
    pts = []
    for g in p.get_drawings():
        r, f = g["rect"], g.get("fill")
        if f and r.y1 < 232 and x0 + 30 < (r.x0 + r.x1) / 2 < x1 and 2.5 < r.width < 5 and np.allclose(f, col, atol=0.02):
            c = ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
            if all(np.hypot(c[0] - q[0], c[1] - q[1]) > 2 for q in pts):
                pts.append(c)
    pts.sort()
    out[name] = [[round(float(10 ** ((x - kx[1]) / kx[0])), 1), round(float(10 ** ((y - ky[1]) / ky[0])), 1)] for x, y in pts]
json.dump(dict(source="Li et al., Nat. Commun. 13, 4901 (2022), Fig. 1a-b; digitised from vector PDF marker centres (OSTI 2470902). "
                      "Axis calibration from tick-label centres; expected digitisation error ~2%.",
               units="T in K, G in MW m^-2 K^-1", data=out), sys.stdout, indent=1)
