#!/usr/bin/env python
"""
T1.1 v2 — audit-corrected analysis: common-index conditional band power + universal-noise envelope.

Fixes over v1 (RESULTS doc, "Revision 2"):
 1. UNIVERSAL INDEX ENFORCED: condition each pulsar's (log10A, gamma) chain on gamma = gamma_c
    (importance reweighting), scanning gamma_c. With gamma common, cross-pulsar comparison is
    reference-frequency independent -- removes the f_yr-vs-low-f ambiguity that flipped v1.
 2. BAND POWER: E = sum_k rho^2 over array bins k=1..5 (T = 16.03 yr). Theory: E ∝ d, log-log
    slope +1 (amplitude slope +0.5).
 3. CENSORED PULSARS INCLUDED: a universal process must appear in EVERY pulsar. Quiet gold
    pulsars set a ceiling K95 = min_i E95_i/d_i on dE/dd. Regression on detections alone cannot
    confirm universality (survivor bias); the envelope can refute or bound it.
 4. CONFOUNDERS: DM covariate added (chromatic leakage grows with DM ~ d); nu-dot Shklovskii-
    corrected (F1_obs contains nu*mu^2*d/c -> distance leaked into v1's covariate).
 5. CALIBRATION: permutation test (distance shuffles) per gamma_c + Sidak over the grid.
 6. VALIDITY: per-pulsar valid band k >= kmin(T_span); J0437 (4.7 yr) excluded; suspect
    parallaxes flagged via implied n_e = DM/d.
"""
import os, glob, re, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0")
ND, PDIR = os.path.join(BASE, "narrowband/noise"), os.path.join(BASE, "narrowband/par")
HD = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/freespectra/ceffyl_data/30f_fs{hd}_100fDMGP_ceffyl")
OUT = os.path.join(ROOT, "data/T1.1_pulsar_timing/derived")
FIG = os.path.join(ROOT, "docs/t1_1_v2_envelope.png")

FYR = 1.0 / (365.25 * 86400)
C_MS = 299792458.0
KPC_M = 3.0856775814913673e19
MASYR_RAD_S = 4.84813681109536e-9 / (365.25 * 86400)
KMAX = 5                      # array band bins k=1..5 (2.0-9.9 nHz), where fs{hd} is measured
GAMMA_GRID = [2.0, 2.5, 3.0, 3.5, 4.0, 13/3, 4.5, 5.0]
H_KERN = 0.15                 # gamma-conditioning kernel width
ESS_MIN = 50
NCOND, NMC, NPERM = 600, 400, 4000
rng = np.random.default_rng(20260724)

# ---------------- free-spectrum HD posterior (per-bin inverse-CDF samplers) ----------------
freqs = np.load(HD + "/freqs.npy")
grid = np.load(HD + "/log10rhogrid.npy")
logden = np.load(HD + "/density.npy")[0]
DF = freqs[0]                                        # = 1/T_array
T_ARR = 1 / DF / (365.25 * 86400)
cdf = np.cumsum(np.exp(logden - logden.max(1, keepdims=True)), 1); cdf /= cdf[:, -1:]
hd_med = np.array([np.interp(0.5, cdf[k], grid) for k in range(KMAX)])
E_GWB = float(np.sum(10 ** (2 * hd_med)))            # HD band power k=1..5 [s^2]
def hd_power_draws(n):
    u = rng.random((n, KMAX))
    return np.sum(10 ** (2 * np.array([np.interp(u[:, k], cdf[k], grid) for k in range(KMAX)]).T), axis=1)

def band_shape(gc, kmin=1):                          # sum f_k^-gc over valid bins [s^-... units folded]
    return np.sum(freqs[kmin - 1:KMAX] ** (-gc))
def phi_base(gc, kmin=1):                            # E = phi_base * A^2
    return (1 / (12 * np.pi ** 2)) * FYR ** (gc - 3) * band_shape(gc, kmin) * DF

# ---------------- load pulsars: chains + par metadata ----------------
PRIMARY = re.compile(r"^[BJ]\d{4}[+-]\d{2,4}$")
def read_par(path):
    d = {}
    for line in open(path):
        t = line.split()
        if t and t[0] in ("PX", "F0", "F1", "DM", "PMELONG", "PMELAT", "PMRA", "PMDEC", "START", "FINISH"):
            try:
                d[t[0]] = float(t[1])
                if t[0] == "PX" and len(t) > 3: d["PX_err"] = float(t[3])
            except ValueError: pass
    return d

parmap = {}
for f in glob.glob(PDIR + "/*.par"):
    parmap.setdefault(os.path.basename(f).split("_")[0], f)

pulsars = {}
for pf in sorted(glob.glob(ND + "/*.pars.txt")):
    psr = os.path.basename(pf).split(".")[0]
    if not PRIMARY.match(psr) or psr not in parmap: continue
    pars = open(pf).read().split()
    iA = [i for i, p in enumerate(pars) if p.endswith("red_noise_log10_A")]
    ig = [i for i, p in enumerate(pars) if p.endswith("red_noise_gamma")]
    if not iA: continue
    iA, ig = iA[0], ig[0]; order = sorted([iA, ig])
    ch = pd.read_csv(ND + f"/{psr}.nb.chain_1.txt", sep=r"\s+", header=None, usecols=order).values
    cm = {order[k]: ch[:, k] for k in range(2)}
    burn = int(0.25 * len(ch))
    par = read_par(parmap[psr])
    if "PX" not in par or "PX_err" not in par or par["PX"] <= 0: continue
    span = (par.get("FINISH", 0) - par.get("START", 0)) / 365.25
    mu = np.hypot(par.get("PMELONG", par.get("PMRA", 0.0)), par.get("PMELAT", par.get("PMDEC", 0.0)))
    d_kpc = 1.0 / par["PX"]
    # Shklovskii-corrected spin-down:  nudot_int = nudot_obs + nu*mu^2*d/c   (obs is negative)
    f1_shk = par["F0"] * (mu * MASYR_RAD_S) ** 2 * d_kpc * KPC_M / C_MS
    f1_int = par["F1"] + f1_shk
    ne = par.get("DM", np.nan) / (d_kpc * 1000)      # implied mean electron density [cm^-3]
    pulsars[psr] = dict(logA=cm[iA][burn:], gam=cm[ig][burn:], PX=par["PX"], PXe=par["PX_err"],
                        d=d_kpc, span=span, DM=par.get("DM", np.nan), F0=par["F0"],
                        F1=par["F1"], F1_int=f1_int, ne=ne,
                        kmin=max(1, int(np.ceil(T_ARR / max(span, 0.5)))))

gold = {p: v for p, v in pulsars.items() if v["PX"] / v["PXe"] > 3}
for p, v in gold.items():
    v["px_suspect"] = bool(v["ne"] > 0.10)           # implied n_e too high -> parallax suspect
    v["band_ok"] = v["kmin"] <= 3                    # span covers at least bins k=3..5
print(f"pulsars={len(pulsars)}  gold={len(gold)}  "
      f"band_ok={sum(v['band_ok'] for v in gold.values())}  "
      f"px_suspect={[p for p,v in gold.items() if v['px_suspect']]}")
print(f"excluded short-span: {[(p, round(v['span'],1)) for p,v in gold.items() if not v['band_ok']]}")

# ---------------- gamma conditioning ----------------
def conditional_A(v, gc, n=NCOND):
    """Samples of log10A conditioned on gamma=gc, + effective sample size."""
    w = np.exp(-0.5 * ((v["gam"] - gc) / H_KERN) ** 2)
    sw = w.sum()
    if sw <= 0: return None, 0.0
    ess = sw ** 2 / np.sum(w ** 2)
    idx = rng.choice(len(w), size=n, p=w / sw)
    return v["logA"][idx], ess

def wquant(x, w, q):
    s = np.argsort(x); x, w = x[s], w[s]
    c = np.cumsum(w); c /= c[-1]
    return np.interp(q, c, x)

# ---------------- per-gamma_c regression with permutation calibration ----------------
results = {"T_array_yr": round(T_ARR, 2), "E_GWB_band": E_GWB, "gamma_grid": {}}
theil_by_gc, det_sets = {}, {}
for gc in GAMMA_GRID:
    pb = phi_base(gc)
    rows = []
    for p, v in gold.items():
        if not v["band_ok"]: continue
        a, ess = conditional_A(v, gc)
        if a is None or ess < ESS_MIN: continue
        # E over k=1..5 under the common-gamma power-law model (the model itself extrapolates
        # below a pulsar's own band; band_ok/kmin gate how far that extrapolation is trusted)
        E = pb * 10 ** (2 * a)
        E_exc_med = pb * 10 ** (2 * np.median(a)) - E_GWB
        p_exc = float(np.mean(E - hd_power_draws(len(E)) > 0))
        rows.append(dict(psr=p, d=v["d"], DM=v["DM"], F0=v["F0"], F1i=v["F1_int"],
                         ess=ess, p_exc=p_exc, E_med=E_exc_med, suspect=v["px_suspect"]))
    R = pd.DataFrame(rows)
    det = R[(R.p_exc > 0.5) & (R.E_med > 0) & (~R.suspect)].copy()
    det_sets[gc] = det
    entry = {"n_compatible": len(R), "n_detected": len(det)}
    if len(det) >= 5:
        x, y = np.log10(det.d.values), np.log10(det.E_med.values)
        ts = stats.theilslopes(y, x)
        sp = stats.spearmanr(x, y)
        # permutation: shuffle distances among the detected set
        perm = np.array([stats.theilslopes(y, x[rng.permutation(len(x))])[0] for _ in range(NPERM)])
        p_perm = float(np.mean(perm >= ts[0]))        # one-sided: theory predicts positive
        # MC slope band over conditional draws + HD draw + parallax noise
        mcs = []
        cond = {r.psr: conditional_A(gold[r.psr], gc)[0] for r in det.itertuples()}
        for _ in range(NMC):
            hd_p = hd_power_draws(1)[0]
            xs, ys = [], []
            for r in det.itertuples():
                a = cond[r.psr][rng.integers(NCOND)]
                E = pb * 10 ** (2 * a) - hd_p
                px = rng.normal(gold[r.psr]["PX"], gold[r.psr]["PXe"])
                if E > 0 and px > 0:
                    xs.append(np.log10(1 / px)); ys.append(np.log10(E))
            if len(xs) >= 5: mcs.append(stats.theilslopes(np.asarray(ys), np.asarray(xs))[0])
        mcs = np.asarray(mcs)
        # covariate OLS: log d + log F0 + log|F1_int| + log DM  (collinearity reported)
        Xc = np.column_stack([np.ones(len(x)), x, np.log10(det.F0), np.log10(np.abs(det.F1i)), np.log10(det.DM)])
        beta = np.linalg.lstsq(Xc, y, rcond=None)[0]
        cc = np.corrcoef(np.column_stack([x, np.log10(det.DM), np.log10(np.abs(det.F1i))]).T)
        entry.update(theil_slope=float(ts[0]), theil_lo=float(ts[2]), theil_hi=float(ts[3]),
                     spearman_rho=float(sp.correlation), spearman_p=float(sp.pvalue),
                     p_perm_onesided=p_perm,
                     mc_slope_med=float(np.median(mcs)), mc_lo68=float(np.percentile(mcs, 16)),
                     mc_hi68=float(np.percentile(mcs, 84)), p_pos=float(np.mean(mcs > 0)),
                     cov_slope_d=float(beta[1]), corr_logd_logDM=float(cc[0, 1]),
                     corr_logd_logF1=float(cc[0, 2]))
        theil_by_gc[gc] = (ts[0], p_perm)
    results["gamma_grid"][f"{gc:.3f}"] = entry

print("\n========== common-gamma conditional regression (theory: slope=+1 in power) ==========")
print(f"{'g_c':>6}{'n_cmp':>6}{'n_det':>6}{'Theil':>8}{'MC med':>8}{'P>0':>6}{'perm p':>8}{'cov_d':>7}{'r(d,DM)':>8}")
for gc in GAMMA_GRID:
    e = results["gamma_grid"][f"{gc:.3f}"]
    if "theil_slope" in e:
        print(f"{gc:6.2f}{e['n_compatible']:6d}{e['n_detected']:6d}{e['theil_slope']:+8.2f}"
              f"{e['mc_slope_med']:+8.2f}{e['p_pos']:6.2f}{e['p_perm_onesided']:8.3f}"
              f"{e['cov_slope_d']:+7.2f}{e['corr_logd_logDM']:8.2f}")
    else:
        print(f"{gc:6.2f}{e['n_compatible']:6d}{e['n_detected']:6d}   (n<5: no fit)")
if theil_by_gc:
    pmin = min(p for _, p in theil_by_gc.values())
    p_glob = 1 - (1 - pmin) ** len(theil_by_gc)
    results["p_min_perm"], results["p_global_sidak"] = float(pmin), float(p_glob)
    print(f"min perm p over grid = {pmin:.3f}  ->  Sidak-corrected global p = {p_glob:.3f}")

# ---------------- universal-noise envelope & ceiling ----------------
GC_ENV = 13/3
pb_env = phi_base(GC_ENV)
env = []
for p, v in gold.items():
    if not v["band_ok"]: continue
    a, ess = conditional_A(v, GC_ENV, n=4000)
    if a is None or ess < 30: continue               # for limits gamma is rail -> ess large anyway
    w_uni = 10 ** a                                  # reweight log-uniform -> uniform-A prior (conservative)
    logA95_lu = np.quantile(a, 0.95)
    logA95_ua = wquant(a, w_uni, 0.95)
    E95 = pb_env * 10 ** (2 * logA95_ua)
    Emed = pb_env * 10 ** (2 * np.median(a))
    p_exc = float(np.mean(pb_env * 10 ** (2 * a[:2000]) - hd_power_draws(2000) > 0))
    env.append(dict(psr=p, d=v["d"], span=v["span"], kmin=v["kmin"], suspect=v["px_suspect"],
                    detected=(p_exc > 0.5) and (Emed > E_GWB), E_med=Emed,
                    E95_uA=E95, E95_lu=pb_env * 10 ** (2 * logA95_lu), K=E95 / v["d"]))
EV = pd.DataFrame(env).sort_values("d").reset_index(drop=True)
quiet = EV[(~EV.detected) & (~EV.suspect)]
kmin_row = quiet.loc[quiet.K.idxmin()]
K95 = float(kmin_row.K)
A_eq_1kpc = 0.5 * np.log10(K95 * 1.0 / pb_env)      # log10 amplitude of universal comp at 1 kpc
loud = EV[EV.detected]
K_loud = float((loud.E_med / loud.d).median()) if len(loud) else np.nan
floor_sp = stats.spearmanr(np.log10(quiet.d), np.log10(quiet.E95_uA))
results["envelope"] = dict(
    gamma_c=GC_ENV, n_envelope=len(EV), n_quiet=len(quiet), n_detected=int(EV.detected.sum()),
    K95_s2_per_kpc=K95, K95_pulsar=str(kmin_row.psr), K95_d_kpc=float(kmin_row.d),
    log10A_universal_max_1kpc=float(A_eq_1kpc), E_GWB_band=E_GWB,
    K_loud_median=K_loud, exclusion_dex=float(np.log10(K_loud / K95)) if K_loud else None,
    floor_spearman_rho=float(floor_sp.correlation), floor_spearman_p=float(floor_sp.pvalue))
print("\n========== universal-noise envelope (gamma_c = 13/3, band k=1-5) ==========")
print(f"HD band power E_GWB = {E_GWB:.2e} s^2  (sqrt = {np.sqrt(E_GWB)*1e9:.0f} ns)")
print(f"quiet gold pulsars (limits): {len(quiet)}  | detected: {int(EV.detected.sum())}")
print(quiet[["psr", "d", "span", "kmin", "E95_uA", "K"]].to_string(index=False,
      float_format=lambda v: f"{v:.3g}"))
print(f"\nCEILING: K95 = {K95:.2e} s^2/kpc  (set by {kmin_row.psr} at {kmin_row.d:.2f} kpc)")
print(f"  -> universal component at 1 kpc: log10 A <= {A_eq_1kpc:.2f}")
print(f"  vs A_gwb = -14.62; vs loud-pulsar K_med = {K_loud:.2e}  "
      f"(exclusion: {np.log10(K_loud/K95):.1f} dex)")
print(f"floor trend (quiet E95 vs d): Spearman rho={floor_sp.correlation:+.2f} (p={floor_sp.pvalue:.2f})"
      f"   [universal process needs floor RISING with d]")

# ---------------- gamma-consistency among native-fit excess pulsars ----------------
det0 = det_sets.get(3.0, pd.DataFrame())
exc_native = [p for p, v in gold.items()
              if v["band_ok"] and phi_base(3.0) * 10 ** (2 * np.median(v["logA"])) > E_GWB
              and np.median(v["logA"]) > -15.5]
if len(exc_native) >= 3:
    gg = np.linspace(0.1, 6.9, 341)
    joint = np.zeros_like(gg); incompat = []
    for p in exc_native:
        kde = stats.gaussian_kde(gold[p]["gam"])
        joint += np.log(kde(gg) + 1e-300)
    gstar = float(gg[np.argmax(joint)])
    for p in exc_native:
        lo, hi = np.percentile(gold[p]["gam"], [2.5, 97.5])
        if not (lo <= gstar <= hi): incompat.append(p)
    results["gamma_consistency"] = dict(n=len(exc_native), gamma_star=gstar,
                                        n_incompatible_95=len(incompat), incompatible=incompat)
    print(f"\ngamma-consistency (n={len(exc_native)} excess pulsars): joint gamma* = {gstar:.2f}; "
          f"{len(incompat)} pulsars exclude gamma* at 95%: {incompat}")

json.dump(results, open(os.path.join(OUT, "t1_1_v2_results.json"), "w"), indent=2, default=str)
EV.to_csv(os.path.join(OUT, "t1_1_v2_envelope_table.csv"), index=False)

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))
a0 = ax[0]
ld = loud
a0.errorbar(ld.d, ld.E_med, fmt="o", ms=7, color="#1f77b4", mec="k", mew=0.4, label="detected excess (median)")
a0.scatter(quiet.d, quiet.E95_uA, marker="v", s=70, facecolors="none", edgecolors="crimson",
           linewidths=1.4, label="quiet: 95% upper limit")
sus = EV[EV.suspect]
a0.scatter(sus.d, sus.E95_uA, marker="x", s=50, color="gray", label="suspect parallax (excl.)")
dd = np.logspace(np.log10(0.2), np.log10(4.5), 50)
a0.plot(dd, K95 * dd, "crimson", lw=2, label=f"universal ceiling $K_{{95}}\\cdot d$")
if np.isfinite(K_loud):
    a0.plot(dd, K_loud * dd, "#1f77b4", ls="--", lw=1.5, alpha=0.8,
            label=f"loud-pulsar trend (excluded as universal, {np.log10(K_loud/K95):.1f} dex above)")
a0.axhline(E_GWB, color="k", ls=":", lw=1.2, label="HD (GWB) band power")
for r in EV.itertuples():
    if r.psr in ("B1937+21", "J1911+1347", "J0740+6620", "J0030+0451"):
        yv = r.E_med if r.detected else r.E95_uA
        a0.annotate(r.psr, (r.d, yv), fontsize=7, xytext=(4, 4), textcoords="offset points")
a0.set_xscale("log"); a0.set_yscale("log")
a0.set_xlabel("distance d [kpc]"); a0.set_ylabel(r"achromatic band power $E$ (k=1–5) [s$^2$]")
a0.set_title(r"Universal-noise envelope ($\gamma_c=13/3$): floor is flat, not rising")
a0.legend(fontsize=7.5, loc="upper left"); a0.grid(alpha=0.3, which="both")

a1 = ax[1]
gcs = [gc for gc in GAMMA_GRID if "theil_slope" in results["gamma_grid"][f"{gc:.3f}"]]
med = [results["gamma_grid"][f"{gc:.3f}"]["mc_slope_med"] for gc in gcs]
lo = [results["gamma_grid"][f"{gc:.3f}"]["mc_lo68"] for gc in gcs]
hi = [results["gamma_grid"][f"{gc:.3f}"]["mc_hi68"] for gc in gcs]
pp = [results["gamma_grid"][f"{gc:.3f}"]["p_perm_onesided"] for gc in gcs]
a1.fill_between(gcs, lo, hi, alpha=0.25, color="seagreen", label="MC 68%")
a1.plot(gcs, med, "o-", color="seagreen", label="Theil–Sen slope (MC median)")
for g, m, p in zip(gcs, med, pp):
    a1.annotate(f"p={p:.2f}", (g, m), fontsize=7, xytext=(0, 8), textcoords="offset points", ha="center")
a1.axhline(1.0, color="b", ls="--", lw=1.5, label="theory: power ∝ d (slope +1)")
a1.axhline(0.0, color="k", lw=1.2, label="null (spin noise)")
a1.set_xlabel(r"common spectral index $\gamma_c$")
a1.set_ylabel(r"slope  d$\log E_{\rm exc}$ / d$\log d$")
a1.set_title(f"Distance slope vs assumed universal index (perm p labels; global p={results.get('p_global_sidak', float('nan')):.2f})")
a1.legend(fontsize=8); a1.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(FIG, dpi=130)
print(f"\nsaved: {FIG}\nsaved: {os.path.join(OUT,'t1_1_v2_results.json')}\nsaved: envelope table csv")
