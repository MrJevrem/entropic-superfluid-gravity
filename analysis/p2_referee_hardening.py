#!/usr/bin/env python
"""
Paper II referee-hardening computations:
 (1) conditioning-kernel sensitivity: K95 anchors at h in {0.10, 0.15, 0.25}
 (2) prior sensitivity: uniform-A vs log-uniform 95% limits
 (3) anchor distance-error budget: J1909-3744 (PPTA) vs J1640+2224 (NG)
 => a defensible K95 +/- uncertainty statement.
"""
import os, sys, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from t2_guard import DERIVED, ROOT

FYR = 1 / (365.25 * 86400)
DF = 1 / (16.03 * 365.25 * 86400)
FK = np.arange(1, 6) * DF
GC = 13 / 3
pb = (1 / (12 * np.pi ** 2)) * FYR ** (GC - 3) * np.sum(FK ** -GC) * DF
rng = np.random.default_rng(20260724)

def wq(x, w, q):
    s = np.argsort(x); x, w = x[s], w[s]
    c = np.cumsum(w); c /= c[-1]
    return np.interp(q, c, x)

def E_of(logA, gam):
    return ((10 ** logA) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3)
            * np.sum(FK[None, :] ** (-np.atleast_1d(gam)[:, None]), 1) * DF)

# ---- load the two anchors ----
PP = os.path.join(ROOT, "data/T1.1_pulsar_timing/ppta_dr3/PPTA-DR3/analysis_codes/data/all/chains")
pars = [l.strip() for l in open(PP + "/commonNoise_pl_nocorr_freegam_DE440_pars.txt") if l.strip()]
ch = np.load(PP + "/chain_commonNoise_pl_nocorr_freegam_DE440.npy"); ch = ch[len(ch) // 4:]
col = {p: i for i, p in enumerate(pars)}
crn = E_of(ch[:, col["gw_pl_nocorr_freegam_log10_A"]], ch[:, col["gw_pl_nocorr_freegam_gamma"]])
anchors = {"J1909-3744(PPTA)": dict(logA=ch[:, col["J1909-3744_red_noise_log10_A"]],
                                    gam=ch[:, col["J1909-3744_red_noise_gamma"]],
                                    common=crn, px=0.86, pxe=0.010)}
ND = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/noise")
pr = open(ND + "/J1640+2224.nb.pars.txt").read().split()
iA = [i for i, p in enumerate(pr) if p.endswith("red_noise_log10_A")][0]
ig = [i for i, p in enumerate(pr) if p.endswith("red_noise_gamma")][0]
o = sorted([iA, ig])
cc = pd.read_csv(ND + "/J1640+2224.nb.chain_1.txt", sep=r"\s+", header=None, usecols=o).values
cm = {o[k]: cc[:, k] for k in range(2)}
E_hd = float(pb * 10 ** (2 * -14.62))
anchors["J1640+2224(NG)"] = dict(logA=cm[iA][len(cc) // 4:], gam=cm[ig][len(cc) // 4:],
                                 common=np.full(len(cc) - len(cc) // 4, E_hd),
                                 px=0.586, pxe=0.183)

out = {}
for name, a in anchors.items():
    d = 1 / a["px"]
    row = dict(d_kpc=round(d, 3), d_frac_err=round(a["pxe"] / a["px"], 3))
    for h in (0.10, 0.15, 0.25):
        w = np.exp(-0.5 * ((a["gam"] - GC) / h) ** 2)
        idx = rng.choice(len(w), 4000, p=w / w.sum())
        aa = a["logA"][idx]
        E_tot = pb * 10 ** (2 * aa) + a["common"][rng.integers(len(a["common"]), size=4000)]
        row[f"K95_uA|h={h}"] = float(wq(E_tot, 10 ** aa, 0.95) / d)
        row[f"K95_lu|h={h}"] = float(np.quantile(E_tot, 0.95) / d)
    out[name] = row

print("=== Paper II hardening numbers ===")
for name, r in out.items():
    print(f"\n[{name}]  d = {r['d_kpc']} kpc (frac err {r['d_frac_err']*100:.1f}%)")
    for h in (0.10, 0.15, 0.25):
        print(f"  h={h}:  K95(uniform-A) = {r[f'K95_uA|h={h}']:.3e}   "
              f"K95(log-uniform) = {r[f'K95_lu|h={h}']:.3e}   "
              f"ratio = {r[f'K95_uA|h={h}']/r[f'K95_lu|h={h}']:.2f}")
j = out["J1909-3744(PPTA)"]
spread_h = max(j[f"K95_uA|h={h}"] for h in (0.10, 0.15, 0.25)) / \
           min(j[f"K95_uA|h={h}"] for h in (0.10, 0.15, 0.25))
print(f"\nJ1909 anchor: h-spread factor = {spread_h:.2f}; distance err = {j['d_frac_err']*100:.1f}%")
print(f"J1640 anchor: distance err = {out['J1640+2224(NG)']['d_frac_err']*100:.0f}% "
      f"-> demote from headline; J1909-anchored value is the quotable ceiling")
json.dump(out, open(os.path.join(DERIVED, "p2_hardening.json"), "w"), indent=1)
print("\nsaved derived/p2_hardening.json")