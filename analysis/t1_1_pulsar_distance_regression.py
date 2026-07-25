#!/usr/bin/env python
"""
T1.1 — Pulsar distance regression  (Entropic-Stochastic Gravity program, RESULTS Sec. 7.4)

THEORY SIGNATURE (discriminator #4): spacetime-diffusion path accumulation makes the
per-pulsar *achromatic* excess red-noise POWER scale with pulsar distance, with a COMMON
spectral index and NO Hellings-Downs correlation. It is therefore separable from:
  - the GW background   -> HD-correlated, distance-independent (a common pedestal), and
  - intrinsic spin noise -> uncorrelated, distance-independent (slope 0 null).

METHOD
  1. Per pulsar i, take the achromatic red-noise posterior (log10 A_i, gamma_i) from the
     NANOGrav 15yr single-pulsar noise analysis (SPNA) MCMC chains. These single-pulsar
     fits ABSORB the common process, so A_i = intrinsic + common (blended).
  2. Subtract the HD-correlated common background at the published posterior level. We work
     at the reference frequency f = 1/yr, where the power-law PSD  P(f)=A^2/(12 pi^2)(f/f_yr)^-g
     reduces to P(f_yr) ∝ A^2  INDEPENDENT of gamma. Hence the excess power at f_yr is
        A_exc,i^2 = A_i^2 - A_gwb^2      (defined where A_i > A_gwb).
     A_gwb from NANOGrav 15yr HD posterior (gamma=13/3): A=2.4e-15 (+0.7/-0.6) ->
     log10 A_gwb ~ N(-14.62, 0.11).
  3. Regress log10 A_exc,i on log10 d_i (d = 1/PX from the timing parallax), controlling for
     the Shannon-Cordes intrinsic-spin-noise proxies log10 nu (=log F0) and log10|nu-dot|
     (=log|F1|). Theory: slope(log A_exc vs log d) = +0.5 (power ∝ d); null (spin noise): 0.
  4. Everything propagated "at the posterior level" by Monte Carlo over the SPNA chains, the
     A_gwb posterior, and the parallax posterior.

Outputs: derived table CSV, results JSON, and a 2-panel figure.
"""
import os, glob, re, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ----------------------------------------------------------------------------- config
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0")
ND   = os.path.join(BASE, "narrowband", "noise")
PDIR = os.path.join(BASE, "narrowband", "par")
OUT  = os.path.join(ROOT, "data/T1.1_pulsar_timing/derived")
FIG  = os.path.join(ROOT, "docs", "t1_1_pulsar_distance_regression.png")
os.makedirs(OUT, exist_ok=True)

PRIMARY   = re.compile(r"^[BJ]\d{4}[+-]\d{2,4}$")   # keep full-array fit; drop *ao/*gbt subsets
LOG10_A_GWB, SIG_A_GWB = -14.62, 0.11               # NANOGrav15 HD, gamma=13/3, A=2.4e-15
PX_SNR_MIN = 3.0                                     # gold parallax subsample
RN_FLOOR   = -15.5                                   # log10 A below this ~ prior rail (undetected RN)
GAMMA_MIN  = 2.0                                     # common-index cut: keep steep, achromatic/spin-GW-like
                                                     # red noise; drop shallow-gamma chromatic contaminants
NSAMP, NMC = 2000, 4000
rng = np.random.default_rng(20260724)

# ----------------------------------------------------------------------------- helpers
def read_par(path):
    d = {}
    for line in open(path):
        t = line.split()
        if not t:
            continue
        k = t[0]
        if k in ("PX", "F0", "F1", "DM", "ELAT"):
            try:
                d[k] = float(t[1])
                if k in ("PX", "F0", "F1") and len(t) > 3:
                    d[k + "_err"] = float(t[3])
            except ValueError:
                pass
    return d

parmap = {}
for f in glob.glob(PDIR + "/*.par"):
    parmap.setdefault(os.path.basename(f).split("_")[0], f)

# ----------------------------------------------------------------------------- load posteriors
rows, samples = [], {}
for pf in sorted(glob.glob(ND + "/*.pars.txt")):
    psr = os.path.basename(pf).split(".")[0]
    if not PRIMARY.match(psr):
        continue
    pars = open(pf).read().split()
    iA = [i for i, p in enumerate(pars) if p.endswith("red_noise_log10_A")]
    ig = [i for i, p in enumerate(pars) if p.endswith("red_noise_gamma")]
    if not iA or psr not in parmap:
        continue
    iA, ig = iA[0], ig[0]
    order = sorted([iA, ig])                          # pandas usecols returns ASCENDING file order
    ch = pd.read_csv(ND + f"/{psr}.nb.chain_1.txt", sep=r"\s+", header=None,
                     usecols=order).values
    colmap = {order[k]: ch[:, k] for k in range(len(order))}   # map by true column index
    c_logA = colmap[iA][int(0.25 * len(ch)):]         # drop burn-in
    c_gam = colmap[ig][int(0.25 * len(ch)):]
    idx = np.linspace(0, len(c_logA) - 1, min(NSAMP, len(c_logA))).astype(int)
    logA, gam = c_logA[idx], c_gam[idx]
    par = read_par(parmap[psr])
    if "PX" not in par or "PX_err" not in par:
        continue
    # P(A_i > A_gwb): significance of excess above the HD background, propagating both posteriors
    p_excess = float(np.mean(logA > rng.normal(LOG10_A_GWB, SIG_A_GWB, size=len(logA))))
    rows.append(dict(psr=psr, logA_med=float(np.median(logA)), gamma_med=float(np.median(gam)),
                     p_excess=p_excess, PX=par["PX"], PX_err=par["PX_err"], F0=par.get("F0", np.nan),
                     F1=par.get("F1", np.nan), DM=par.get("DM", np.nan)))
    samples[psr] = (logA, gam)

df = pd.DataFrame(rows)
df["PX_snr"] = df.PX / df.PX_err
df["dist_kpc"] = 1.0 / df.PX
# per-pulsar median excess amplitude (for plotting / table); NaN if median below GWB
A_gwb2_med = (10 ** LOG10_A_GWB) ** 2
def med_excess(psr):
    # subtract at the median amplitude (no per-sample conditioning -> unbiased); NaN if median below GWB
    logA, _ = samples[psr]
    exc = (10 ** np.median(logA)) ** 2 - A_gwb2_med
    return float(np.sqrt(exc)) if exc > 0 else np.nan
df["A_exc_med"] = df.psr.map(med_excess)
df["logA_exc_med"] = np.log10(df.A_exc_med)
# excess: red-noise amplitude posterior more likely than not above the HD background
df["excess_detected"] = (df.p_excess > 0.5) & (df.logA_med > RN_FLOOR)

df = df.sort_values("PX_snr", ascending=False).reset_index(drop=True)
df.to_csv(os.path.join(OUT, "t1_1_pulsar_table.csv"), index=False)

gold = df[df.PX_snr > PX_SNR_MIN].copy()
detected = gold[gold.excess_detected].copy()
print(f"pulsars loaded (primary): {len(df)}")
print(f"gold parallax (PX/sig>{PX_SNR_MIN}): {len(gold)}")
print(f"gold with detected excess above HD background: {len(detected)}")
print("\nDetected-excess sample (theory sample):")
print(detected[["psr", "PX_snr", "dist_kpc", "logA_med", "gamma_med", "p_excess", "logA_exc_med"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ----------------------------------------------------------------------------- MC regression
def mc_regression(pool, use_cov):
    """Posterior-level MC. Returns slope array for log10 A_exc vs log10 d."""
    psrs = pool.psr.values
    PX, PXe = pool.PX.values, pool.PX_err.values
    lF0 = np.log10(pool.F0.values)
    lF1 = np.log10(np.abs(pool.F1.values))
    slopes = []
    for _ in range(NMC):
        A_gwb2 = (10 ** rng.normal(LOG10_A_GWB, SIG_A_GWB)) ** 2
        xs, ys, c0, c1 = [], [], [], []
        for k, psr in enumerate(psrs):
            logA, _ = samples[psr]
            j = rng.integers(len(logA))
            exc = (10 ** logA[j]) ** 2 - A_gwb2
            if exc <= 0:
                continue
            px = rng.normal(PX[k], PXe[k])
            if px <= 0:
                continue
            xs.append(np.log10(1.0 / px))
            ys.append(0.5 * np.log10(exc))            # log10 A_exc = 0.5 log10(power)
            c0.append(lF0[k]); c1.append(lF1[k])
        if len(xs) < 5:
            continue
        xs, ys = np.asarray(xs), np.asarray(ys)
        if use_cov:
            X = np.column_stack([np.ones_like(xs), xs, np.asarray(c0), np.asarray(c1)])
            beta = np.linalg.lstsq(X, ys, rcond=None)[0]
            slopes.append(beta[1])
        else:
            slopes.append(np.polyfit(xs, ys, 1)[0])
    return np.asarray(slopes)

def summ(a):
    return dict(median=float(np.median(a)), lo68=float(np.percentile(a, 16)),
                hi68=float(np.percentile(a, 84)), lo95=float(np.percentile(a, 2.5)),
                hi95=float(np.percentile(a, 97.5)), p_pos=float(np.mean(a > 0)), n=int(len(a)))

slope_biv = mc_regression(detected, use_cov=False)
slope_cov = mc_regression(detected, use_cov=True)

# COMMON-INDEX sample: the theory requires a common (steep) spectral index, so drop shallow-gamma
# (likely chromatic / DM / system) contaminants. This is the theory-faithful sample.
common = detected[detected.gamma_med > GAMMA_MIN].copy()
slope_ci_biv = mc_regression(common, use_cov=False)
slope_ci_cov = mc_regression(common, use_cov=True)

# point-estimate OLS (headline p-value) on median values
m = detected.dropna(subset=["logA_exc_med"])
ols = stats.linregress(np.log10(m.dist_kpc), m.logA_exc_med)
mc_ = common.dropna(subset=["logA_exc_med"])
ols_c = stats.linregress(np.log10(mc_.dist_kpc), mc_.logA_exc_med)

# selection-bias check: linear excess-POWER vs distance on the FULL gold sample (incl. below-GWB)
def mc_linear_power(pool):
    psrs = pool.psr.values; PX, PXe = pool.PX.values, pool.PX_err.values
    sl = []
    for _ in range(NMC):
        A_gwb2 = (10 ** rng.normal(LOG10_A_GWB, SIG_A_GWB)) ** 2
        xs, ys = [], []
        for k, psr in enumerate(psrs):
            logA, _ = samples[psr]
            j = rng.integers(len(logA))
            px = rng.normal(PX[k], PXe[k])
            if px <= 0:
                continue
            xs.append(1.0 / px)
            ys.append(((10 ** logA[j]) ** 2 - A_gwb2) / 1e-27)  # excess power, scaled units
        if len(xs) < 5:
            continue
        sl.append(np.polyfit(xs, ys, 1)[0])
    return np.asarray(sl)

slope_lin = mc_linear_power(gold)
slope_ci_lin = mc_linear_power(common)   # robust primary: excess POWER vs distance, common-index sample

# ---- robustness diagnostics on the common-index sample (median points) ----
xdat = np.log10(mc_.dist_kpc).values
ydat = mc_.logA_exc_med.values
sp = stats.spearmanr(xdat, ydat)
kt = stats.kendalltau(xdat, ydat)
jk = np.array([np.polyfit(np.delete(xdat, i), np.delete(ydat, i), 1)[0] for i in range(len(xdat))])
i_infl = int(np.argmax(np.abs(jk - np.polyfit(xdat, ydat, 1)[0])))
psr_infl = mc_.psr.values[i_infl]
# selection probe: is "excess detected" correlated with distance across the gold sample?
sel_sp = stats.spearmanr(np.log10(gold.dist_kpc), gold.excess_detected.astype(float))
diagnostics = dict(
    spearman=dict(rho=float(sp.correlation), p=float(sp.pvalue)),
    kendall=dict(tau=float(kt.correlation), p=float(kt.pvalue)),
    jackknife_slope=dict(min=float(jk.min()), max=float(jk.max()),
                         most_influential=psr_infl, slope_without_it=float(jk[i_infl])),
    selection_probe=dict(spearman_rho_dist_vs_detected=float(sel_sp.correlation),
                         p=float(sel_sp.pvalue)),
)

results = dict(
    n_primary=len(df), n_gold=len(gold), n_detected=len(detected), n_common_index=len(common),
    A_gwb=dict(log10=LOG10_A_GWB, sigma=SIG_A_GWB, gamma="13/3", A="2.4e-15"),
    slope_bivariate=summ(slope_biv),
    slope_with_covariates=summ(slope_cov),
    slope_commonidx_bivariate=summ(slope_ci_biv),
    slope_commonidx_covariates=summ(slope_ci_cov),
    slope_linear_power_fullgold=summ(slope_lin),
    slope_linear_power_commonidx=summ(slope_ci_lin),
    ols_point_estimate=dict(slope=float(ols.slope), stderr=float(ols.stderr),
                            intercept=float(ols.intercept), rvalue=float(ols.rvalue),
                            pvalue_twosided=float(ols.pvalue)),
    ols_common_index=dict(slope=float(ols_c.slope), stderr=float(ols_c.stderr),
                          rvalue=float(ols_c.rvalue), pvalue_twosided=float(ols_c.pvalue)),
    gamma_excess=dict(median=float(detected.gamma_med.median()),
                      std=float(detected.gamma_med.std()),
                      values={r.psr: round(r.gamma_med, 2) for r in detected.itertuples()}),
    theory_expected_slope=0.5, null_slope=0.0,
    diagnostics=diagnostics,
)
json.dump(results, open(os.path.join(OUT, "t1_1_results.json"), "w"), indent=2)

print("\n================= RESULTS =================")
print(f"OLS (median pts):  slope = {ols.slope:+.3f} +/- {ols.stderr:.3f}  "
      f"(r={ols.rvalue:+.2f}, two-sided p={ols.pvalue:.3f})")
print(f"MC bivariate:      slope = {results['slope_bivariate']['median']:+.3f}  "
      f"[68% {results['slope_bivariate']['lo68']:+.3f},{results['slope_bivariate']['hi68']:+.3f}]  "
      f"P(>0)={results['slope_bivariate']['p_pos']:.3f}")
print(f"MC +covariates:    slope = {results['slope_with_covariates']['median']:+.3f}  "
      f"[68% {results['slope_with_covariates']['lo68']:+.3f},{results['slope_with_covariates']['hi68']:+.3f}]  "
      f"P(>0)={results['slope_with_covariates']['p_pos']:.3f}")
print(f"-- common-index sample (gamma>{GAMMA_MIN}, n={len(common)}) --")
print(f"OLS (median pts):  slope = {ols_c.slope:+.3f} +/- {ols_c.stderr:.3f}  "
      f"(r={ols_c.rvalue:+.2f}, two-sided p={ols_c.pvalue:.3f})")
print(f"MC bivariate:      slope = {results['slope_commonidx_bivariate']['median']:+.3f}  "
      f"[68% {results['slope_commonidx_bivariate']['lo68']:+.3f},{results['slope_commonidx_bivariate']['hi68']:+.3f}]  "
      f"P(>0)={results['slope_commonidx_bivariate']['p_pos']:.3f}")
print(f"MC +covariates:    slope = {results['slope_commonidx_covariates']['median']:+.3f}  "
      f"[68% {results['slope_commonidx_covariates']['lo68']:+.3f},{results['slope_commonidx_covariates']['hi68']:+.3f}]  "
      f"P(>0)={results['slope_commonidx_covariates']['p_pos']:.3f}")
print(f"MC linear power (full gold, contaminated):    P(slope>0)="
      f"{results['slope_linear_power_fullgold']['p_pos']:.3f}")
print(f"MC linear power (common-index, ROBUST PRIMARY): P(slope>0)="
      f"{results['slope_linear_power_commonidx']['p_pos']:.3f}")
print(f"excess-sample gamma: median={results['gamma_excess']['median']:.2f} "
      f"std={results['gamma_excess']['std']:.2f}  (theory: common index)")
print("-- robustness (common-index median pts) --")
print(f"Spearman rho={sp.correlation:+.3f} (p={sp.pvalue:.3f}) | Kendall tau={kt.correlation:+.3f} (p={kt.pvalue:.3f})")
print(f"jackknife slope range [{jk.min():+.3f},{jk.max():+.3f}]; most influential={psr_infl} -> slope {jk[i_infl]:+.3f}")
print(f"selection probe (dist vs detected, gold): Spearman rho={sel_sp.correlation:+.3f} (p={sel_sp.pvalue:.3f})")
print("Theory slope=+0.5 | Null (spin noise) slope=0")

# ----------------------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
mm = detected.dropna(subset=["logA_exc_med"])
xe = mm.PX_err / (mm.PX * np.log(10))                       # d error in dex
sc = ax[0].scatter(np.log10(mm.dist_kpc), mm.logA_exc_med, c=mm.gamma_med,
                   cmap="viridis", s=60, zorder=3, edgecolor="k", linewidth=0.4)
ax[0].errorbar(np.log10(mm.dist_kpc), mm.logA_exc_med, xerr=xe, fmt="none",
               ecolor="gray", alpha=0.5, zorder=2)
xx = np.linspace(np.log10(mm.dist_kpc).min(), np.log10(mm.dist_kpc).max(), 50)
ax[0].plot(xx, ols.intercept + ols.slope * xx, "r-", lw=2,
           label=f"all-detected fit={ols.slope:+.2f}")
ax[0].plot(xx, ols_c.intercept + ols_c.slope * xx, "g-", lw=2,
           label=f"common-$\\gamma$ fit={ols_c.slope:+.2f}")
low = mm[mm.gamma_med <= GAMMA_MIN]
ax[0].scatter(np.log10(low.dist_kpc), low.logA_exc_med, s=150, facecolors="none",
              edgecolors="red", linewidths=1.4, zorder=4, label=f"$\\gamma<{GAMMA_MIN}$ (chromatic?)")
xc, yc = np.log10(mc_.dist_kpc).mean(), mc_.logA_exc_med.mean()
ax[0].plot(xx, yc + 0.5 * (xx - xc), "b--", lw=1.2, alpha=0.7, label="theory slope=+0.5")
for r in mm.itertuples():
    ax[0].annotate(r.psr, (np.log10(r.dist_kpc), r.logA_exc_med), fontsize=6,
                   xytext=(3, 3), textcoords="offset points")
fig.colorbar(sc, ax=ax[0], label=r"red-noise index $\gamma$")
ax[0].set_xlabel(r"$\log_{10}(d\,/\,\mathrm{kpc})$")
ax[0].set_ylabel(r"$\log_{10} A_{\rm exc}$  (achromatic excess, at $f=1/\rm yr$)")
ax[0].set_title("T1.1 pulsar distance regression (NANOGrav 15yr)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

ax[1].hist(slope_cov, bins=50, alpha=0.55, label="all-detected +cov", color="orange", density=True)
ax[1].hist(slope_ci_cov, bins=50, alpha=0.55, label=r"common-$\gamma$ +cov", color="seagreen", density=True)
ax[1].axvline(0, color="k", lw=1.5, label="null (spin noise)")
ax[1].axvline(0.5, color="b", ls="--", lw=1.5, label="theory (+0.5)")
ax[1].set_xlabel("regression slope  d(log $A_{\\rm exc}$)/d(log $d$)")
ax[1].set_ylabel("posterior density")
ax[1].set_title(f"slope posterior  (common-$\\gamma$ P>0 = {results['slope_commonidx_covariates']['p_pos']:.2f})")
ax[1].legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIG, dpi=130)
print(f"\nsaved: {FIG}")
print(f"saved: {os.path.join(OUT,'t1_1_pulsar_table.csv')}")
print(f"saved: {os.path.join(OUT,'t1_1_results.json')}")
