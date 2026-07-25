#!/usr/bin/env python
"""
T1.1 robustness: model-independent HD subtraction via the free-spectrum posterior.

The primary analysis (t1_1_pulsar_distance_regression.py) subtracts the HD background as a
single gamma=13/3 power-law amplitude at f=1/yr. This script instead subtracts the HD
common process PER FREQUENCY using the NANOGrav 15yr HD free-spectrum posterior
(ceffyl 30f_fs{hd}_100fDMGP), which assumes no spectral shape.

Convention (verified): enterprise rho_k^2 = (A^2/12pi^2) f_yr^(g-3) f_k^(-g) * df   [df = 1/T],
and the free spectrum stores rho_hd,k^2 = 10^(2 log10rho). ceffyl 'density' is the LOG pdf.

Per-pulsar excess: E_i = sum over the GWB-dominated band of (rho_i,k^2 - rho_hd,k^2).
Effective excess amplitude A_exc = sqrt(E_i). Regress log10 A_exc on log10 d, common-index
sample (same pulsars as the primary), MC over pulsar chains + HD KDE + parallax.
"""
import os, glob, re, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0")
ND   = os.path.join(BASE, "narrowband", "noise")
HD   = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/freespectra/ceffyl_data/30f_fs{hd}_100fDMGP_ceffyl")
OUT  = os.path.join(ROOT, "data/T1.1_pulsar_timing/derived")
TABLE = os.path.join(OUT, "t1_1_pulsar_table.csv")
PRIM  = os.path.join(OUT, "t1_1_results.json")

FYR = 1.0 / (365.25 * 86400)
NSAMP, NMC = 2000, 4000
BANDS = {"k1-5": 5, "k1-8": 8, "k1-14": 14}     # GWB-dominated bands to test
rng = np.random.default_rng(20260724)

# ---- HD free-spectrum posterior: per-bin inverse-CDF samplers ----
freqs = np.load(HD + "/freqs.npy")
grid  = np.load(HD + "/log10rhogrid.npy")
logden = np.load(HD + "/density.npy")[0]                 # (30, Ngrid) log-pdf
df_spacing = freqs[0]                                    # uniform: df = f_1 = 1/T
cdf = np.cumsum(np.exp(logden - logden.max(axis=1, keepdims=True)), axis=1)
cdf /= cdf[:, -1:]                                       # (30, Ngrid) normalized CDF per bin
def sample_hd(nbins):
    """Draw one HD spectrum: log10 rho_hd for bins 0..nbins-1."""
    u = rng.random(nbins)
    return np.array([np.interp(u[k], cdf[k], grid) for k in range(nbins)])

# ---- pulsar red-noise posteriors (same corrected loader as the primary) ----
PRIMARY = re.compile(r"^[BJ]\d{4}[+-]\d{2,4}$")
samples = {}
for pf in sorted(glob.glob(ND + "/*.pars.txt")):
    psr = os.path.basename(pf).split(".")[0]
    if not PRIMARY.match(psr):
        continue
    pars = open(pf).read().split()
    iA = [i for i, p in enumerate(pars) if p.endswith("red_noise_log10_A")]
    ig = [i for i, p in enumerate(pars) if p.endswith("red_noise_gamma")]
    if not iA:
        continue
    iA, ig = iA[0], ig[0]
    order = sorted([iA, ig])
    ch = pd.read_csv(ND + f"/{psr}.nb.chain_1.txt", sep=r"\s+", header=None, usecols=order).values
    colmap = {order[k]: ch[:, k] for k in range(len(order))}
    cL, cG = colmap[iA][int(0.25 * len(ch)):], colmap[ig][int(0.25 * len(ch)):]
    idx = np.linspace(0, len(cL) - 1, min(NSAMP, len(cL))).astype(int)
    samples[psr] = (cL[idx], cG[idx])

# ---- common-index sample = exactly the primary's theory sample ----
tab = pd.read_csv(TABLE)
sample = tab[(tab.PX_snr > 3) & (tab.excess_detected) & (tab.gamma_med > 2.0)].copy()
sample = sample[sample.psr.isin(samples)].reset_index(drop=True)
print(f"common-index sample (same as primary): n={len(sample)}")
print("  " + ", ".join(sample.psr))

def rho_i2(logA, gam, fk):
    """enterprise per-bin power for a pulsar power-law, at frequencies fk."""
    return (10 ** logA) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3) * fk ** (-gam) * df_spacing

def summ(a):
    a = np.asarray(a)
    return dict(median=float(np.median(a)), lo68=float(np.percentile(a, 16)),
                hi68=float(np.percentile(a, 84)), p_pos=float(np.mean(a > 0)), n=int(len(a)))

def mc_freespectrum(nbins, use_cov):
    fk = freqs[:nbins]
    psrs = sample.psr.values
    PX, PXe = sample.PX.values, sample.PX_err.values
    lF0, lF1 = np.log10(sample.F0.values), np.log10(np.abs(sample.F1.values))
    slopes, ndrop = [], 0
    for _ in range(NMC):
        rho_hd2 = 10 ** (2 * sample_hd(nbins))              # common HD spectrum this draw
        xs, ys, c0, c1 = [], [], [], []
        for k, psr in enumerate(psrs):
            logA, gam = samples[psr]
            j = rng.integers(len(logA))
            E = np.sum(rho_i2(logA[j], gam[j], fk) - rho_hd2)
            if E <= 0:
                continue
            px = rng.normal(PX[k], PXe[k])
            if px <= 0:
                continue
            xs.append(np.log10(1.0 / px)); ys.append(0.5 * np.log10(E))
            c0.append(lF0[k]); c1.append(lF1[k])
        if len(xs) < 5:
            ndrop += 1; continue
        xs, ys = np.asarray(xs), np.asarray(ys)
        if use_cov:
            X = np.column_stack([np.ones_like(xs), xs, np.asarray(c0), np.asarray(c1)])
            slopes.append(np.linalg.lstsq(X, ys, rcond=None)[0][1])
        else:
            slopes.append(np.polyfit(xs, ys, 1)[0])
    return np.asarray(slopes), ndrop

# median-point excess (band k1-5) for a non-parametric check + table
def median_excess_amp(psr, nbins):
    logA, gam = samples[psr]
    rho_hd2 = 10 ** (2 * np.array([np.interp(0.5, cdf[k], grid) for k in range(nbins)]))
    E = np.sum(rho_i2(np.median(logA), np.median(gam), freqs[:nbins]) - rho_hd2)
    return 0.5 * np.log10(E) if E > 0 else np.nan

print("\n================= FREE-SPECTRUM ROBUSTNESS =================")
results = {"method": "per-frequency HD free-spectrum subtraction (30f_fs{hd}_100fDMGP)",
           "n_sample": int(len(sample)), "bands": {}}
for label, nb in BANDS.items():
    sb, sc = mc_freespectrum(nb, use_cov=False), None
    scov, _ = mc_freespectrum(nb, use_cov=True)
    # median-point OLS + rank test on band k1-5-style medians
    sample[f"lAexc_{label}"] = [median_excess_amp(p, nb) for p in sample.psr]
    mm = sample.dropna(subset=[f"lAexc_{label}"])
    ols = stats.linregress(np.log10(mm.dist_kpc), mm[f"lAexc_{label}"])
    spr = stats.spearmanr(np.log10(mm.dist_kpc), mm[f"lAexc_{label}"])
    results["bands"][label] = dict(
        nbins=nb, n_points_median=int(len(mm)),
        ols_slope=float(ols.slope), ols_stderr=float(ols.stderr), ols_p=float(ols.pvalue),
        spearman_rho=float(spr.correlation), spearman_p=float(spr.pvalue),
        mc_bivariate=summ(sb[0]), mc_covariates=summ(scov))
    r = results["bands"][label]
    print(f"\n[{label}] n_med={len(mm)}  OLS slope={ols.slope:+.3f}±{ols.stderr:.3f} (p={ols.pvalue:.3f})  "
          f"Spearman rho={spr.correlation:+.3f} (p={spr.pvalue:.3f})")
    print(f"        MC bivariate slope={r['mc_bivariate']['median']:+.3f} "
          f"[68%{r['mc_bivariate']['lo68']:+.2f},{r['mc_bivariate']['hi68']:+.2f}] P(>0)={r['mc_bivariate']['p_pos']:.3f}")
    print(f"        MC +cov     slope={r['mc_covariates']['median']:+.3f} P(>0)={r['mc_covariates']['p_pos']:.3f}")

# ---- comparison to the power-law primary ----
prim = json.load(open(PRIM))
pc = prim["slope_commonidx_bivariate"]; pcov = prim["slope_commonidx_covariates"]
print("\n----------------- COMPARISON (common-index sample) -----------------")
print(f"{'method':28}{'slope(biv)':>12}{'P(>0)':>8}{'slope(+cov)':>13}{'P(>0)':>8}")
print(f"{'power-law @f_yr (primary)':28}{pc['median']:+12.3f}{pc['p_pos']:>8.3f}"
      f"{pcov['median']:+13.3f}{pcov['p_pos']:>8.3f}")
for label in BANDS:
    b = results["bands"][label]
    print(f"{'free-spec '+label:28}{b['mc_bivariate']['median']:+12.3f}{b['mc_bivariate']['p_pos']:>8.3f}"
          f"{b['mc_covariates']['median']:+13.3f}{b['mc_covariates']['p_pos']:>8.3f}")
results["primary_powerlaw"] = dict(bivariate=pc, covariates=pcov)
json.dump(results, open(os.path.join(OUT, "t1_1_freespectrum_results.json"), "w"), indent=2)
print("\nsaved:", os.path.join(OUT, "t1_1_freespectrum_results.json"))
