#!/usr/bin/env python
"""Comparison figure: power-law-at-f_yr vs model-independent free-spectrum HD subtraction.
Shows the tentative distance trend is not robust to the subtraction method."""
import os, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ND = f"{ROOT}/data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/noise"
HD = f"{ROOT}/data/T1.1_pulsar_timing/nanograv_15yr/freespectra/ceffyl_data/30f_fs{{hd}}_100fDMGP_ceffyl"
FYR = 1 / (365.25 * 86400)
freqs = np.load(HD + "/freqs.npy"); grid = np.load(HD + "/log10rhogrid.npy"); logden = np.load(HD + "/density.npy")[0]
cdf = np.cumsum(np.exp(logden - logden.max(1, keepdims=True)), 1); cdf /= cdf[:, -1:]
df = freqs[0]
rho_hd2 = 10 ** (2 * np.array([np.interp(0.5, cdf[k], grid) for k in range(5)]))

tab = pd.read_csv(f"{ROOT}/data/T1.1_pulsar_timing/derived/t1_1_pulsar_table.csv")
s = tab[(tab.PX_snr > 3) & (tab.excess_detected) & (tab.gamma_med > 2)].copy()

def med_AG(psr):
    pars = open(f"{ND}/{psr}.nb.pars.txt").read().split()
    iA = [i for i, p in enumerate(pars) if p.endswith("red_noise_log10_A")][0]
    ig = [i for i, p in enumerate(pars) if p.endswith("red_noise_gamma")][0]
    o = sorted([iA, ig]); ch = pd.read_csv(f"{ND}/{psr}.nb.chain_1.txt", sep=r"\s+", header=None, usecols=o).values
    cm = {o[k]: ch[:, k] for k in range(2)}; return np.median(cm[iA]), np.median(cm[ig])

fs = []
for r in s.itertuples():
    lA, g = med_AG(r.psr)
    E = np.sum((10 ** lA) ** 2 / (12 * np.pi ** 2) * FYR ** (g - 3) * freqs[:5] ** (-g) * df - rho_hd2)
    fs.append(0.5 * np.log10(E) if E > 0 else np.nan)
s["fs_exc"] = fs
s["logd"] = np.log10(s.dist_kpc)

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for a, (col, ttl, ylab) in zip(ax, [("logA_exc_med", "Power-law @ $f=1/$yr  (primary)", r"$\log_{10}A_{\rm exc}$ at $f_{\rm yr}$"),
                                     ("fs_exc", "Free-spectrum, per-frequency  (model-independent)", r"$\log_{10}\sqrt{E_{\rm exc}}$ (integrated $k$=1-5)")]):
    d = s.dropna(subset=[col])
    ols = stats.linregress(d.logd, d[col]); sp = stats.spearmanr(d.logd, d[col])
    a.scatter(d.logd, d[col], c=d.gamma_med, cmap="viridis", s=70, edgecolor="k", linewidth=0.4, zorder=3)
    xx = np.linspace(d.logd.min(), d.logd.max(), 40)
    a.plot(xx, ols.intercept + ols.slope * xx, "r-", lw=2, label=f"slope={ols.slope:+.2f}, Spearman p={sp.pvalue:.2f}")
    for r in d.itertuples():
        c = "red" if r.psr in ("J0030+0451", "B1937+21") else "0.3"
        a.annotate(r.psr, (r.logd, getattr(r, col)), fontsize=6.5, color=c, xytext=(3, 3), textcoords="offset points")
    a.set_xlabel(r"$\log_{10}(d/{\rm kpc})$"); a.set_ylabel(ylab); a.set_title(ttl)
    a.legend(fontsize=9); a.grid(alpha=0.3)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(s.gamma_med.min(), s.gamma_med.max()))
fig.colorbar(sm, ax=ax, label=r"red-noise index $\gamma$", fraction=0.025)
fig.suptitle("T1.1 robustness: the distance trend does not survive model-independent HD subtraction", fontsize=12)
out = f"{ROOT}/docs/t1_1_subtraction_method_comparison.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print("saved:", out)
print(f"J0030+0451: nearest pulsar (d={s[s.psr=='J0030+0451'].dist_kpc.values[0]:.2f} kpc), "
      f"gamma={s[s.psr=='J0030+0451'].gamma_med.values[0]:.2f} -> flips from low (PL) to high (FS)")
