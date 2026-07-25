#!/usr/bin/env python
"""
T2.2 — sightline-entropy regression (B4's direct test; plan §4).

Question: does ACHROMATIC per-pulsar excess track the thermodynamic entropy of the traversed
medium at fixed DM and distance? B4 predicts yes (frozen: slope=1 in logE-logs, achromatic-
specific); all other branches predict no residual dependence beyond chromatic ISM leakage.

Responses: achromatic E_RN from the JOINT fits (NG curn core, 67 psr; PPTA CURN chain, 30 psr);
chromatic control E_DMGP from the PPTA chain (NG models DM as DMX in the timing fit — no
chromatic GP in its core; deviation logged: chromatic control is PPTA-only).
Covariates: DM (psrcat/par), d (psrcat PX; gold subset primary), EM from the Finkbeiner 2003
Halpha composite (EM = 2.75 T4^0.9 I_Ha), <n_e> = EM/DM, s_path = DM*ln(T_wim^1.5/<n_e>)
(two-phase WIM proxy, frozen), HII indicator = EM above 90th percentile (pygedm build failed
on this platform — NE2001 clump tagging replaced by the EM flag; deviation logged).
Statistics: partial Spearman (controls: log DM, log d), logistic loud/quiet, permutation-
calibrated (1e4), Sidak over the 3 pre-registered proxies; achromatic-vs-chromatic differential.
"""
import os, sys, json, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import h5py, healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats
from t2_guard import enforce, DERIVED, ROOT

HAL = os.path.join(ROOT, "data/T2/lambda_halpha_fwhm06_0512.fits")
freezes, _ = enforce(staged=[HAL])
fz = freezes[-1]
B4P = fz["branches"]["B4"]["t22_prediction"]

FYR = 1 / (365.25 * 86400)
DF = 1 / (16.03 * 365.25 * 86400)
FK = np.arange(1, 6) * DF
def E_of(logA, gam):
    gam = np.atleast_1d(gam)
    return ((10 ** np.atleast_1d(logA)) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3)
            * np.sum(FK[None, :] ** (-gam[:, None]), 1) * DF)

# ---------------- responses: joint achromatic (+ PPTA chromatic control) ----------------
CORE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/gwb_cores/curn_14f_pl_vg.core")
with h5py.File(CORE) as f:
    pars = [p.decode() if isinstance(p, bytes) else str(p) for p in f["params"][:]]
    ch = f["chain"][:][int(f["metadata/burn"][()]):]
col = {p: i for i, p in enumerate(pars)}
NGP = os.path.join(ROOT, "data/T1.1_pulsar_timing/ppta_dr3/PPTA-DR3/analysis_codes/data/all")
pp_pars = [l.strip() for l in open(NGP + "/chains/commonNoise_pl_nocorr_freegam_DE440_pars.txt") if l.strip()]
pch = np.load(NGP + "/chains/chain_commonNoise_pl_nocorr_freegam_DE440.npy")
pch = pch[len(pch) // 4:]
pcol = {p: i for i, p in enumerate(pp_pars)}

v4 = json.load(open(os.path.join(ROOT, "data/T1.1_pulsar_timing/derived/t1_1_v4_results.json")))
loud_ng = set(v4["ng_census"]["detected"]); pdet_ng = v4["ng_census"]["pdet"]
loud_pp = {"J1643-1224", "J1939+2134"}

rows = []
for p in sorted(set(k.split("_red_noise")[0] for k in pars if "_red_noise_log10_A" in k)):
    E = float(np.median(E_of(ch[:, col[f"{p}_red_noise_log10_A"]], ch[:, col[f"{p}_red_noise_gamma"]])))
    rows.append(dict(psr=p, array="NG", logE=np.log10(E), loud=p in loud_ng, logE_chrom=np.nan))
for p in sorted(set(k.split("_red_noise")[0] for k in pp_pars if "_red_noise_log10_A" in k)):
    E = float(np.median(E_of(pch[:, pcol[f"{p}_red_noise_log10_A"]], pch[:, pcol[f"{p}_red_noise_gamma"]])))
    Ec = np.nan
    if f"{p}_dm_gp_log10_A" in pcol:
        Ec = float(np.median(E_of(pch[:, pcol[f"{p}_dm_gp_log10_A"]], pch[:, pcol[f"{p}_dm_gp_gamma"]])))
    rows.append(dict(psr=p, array="PPTA", logE=np.log10(E), loud=p in loud_pp,
                     logE_chrom=np.log10(Ec) if np.isfinite(Ec) else np.nan))
R = pd.DataFrame(rows)

# ---------------- sightline covariates ----------------
ALIAS = {"B1937+21": "J1939+2134", "B1855+09": "J1857+0943"}
def compact_err(valstr, errstr):
    m = re.match(r"^[-+]?\d+(\.(\d+))?([eE]([-+]?\d+))?$", valstr)
    if not m or not re.match(r"^\d+$", errstr): return None
    nd = len(m.group(2) or ""); exp = int(m.group(4) or 0)
    return int(errstr) * (10.0 ** (exp - nd))
cat, cur = {}, {}
for raw in open(os.path.join(ROOT, "data/T1.1_pulsar_timing/atnf_psrcat/psrcat_tar/psrcat.db"),
                encoding="latin-1"):
    l = raw.rstrip("\n")
    if l.startswith("@"):
        if cur:
            for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = dict(cur)
        cur = {}; continue
    if not l.strip() or l.startswith("#"): continue
    t = l.split()
    if t[0] == "PSRJ": cur["_J"] = t[1]
    elif t[0] == "PSRB": cur["_B"] = t[1]
    elif t[0] in ("RAJ", "DECJ"):
        cur[t[0]] = t[1]                                  # sexagesimal string
    elif t[0] in ("DM", "DIST_DM"):
        try: cur[t[0]] = float(t[1])
        except ValueError: pass
    elif t[0] == "PX" and len(t) >= 3:
        try:
            cur["PX"] = float(t[1]); e = compact_err(t[1], t[2])
            if e: cur["PXe"] = e
        except ValueError: pass
if cur:
    for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = dict(cur)

halpha = hp.read_map(HAL, verbose=False) if "verbose" in hp.read_map.__code__.co_varnames \
    else hp.read_map(HAL)
NSIDE = hp.get_nside(halpha)
T4 = 0.8                                                    # T_WIM = 8000 K (frozen)
meta = []
for r in R.itertuples():
    c = cat.get(r.psr) or cat.get(ALIAS.get(r.psr, ""), None)
    if not c or "RAJ" not in c or "DECJ" not in c or "DM" not in c:
        meta.append(dict(DM=np.nan, d=np.nan)); continue
    sc = SkyCoord(ra=c["RAJ"], dec=c["DECJ"], unit=(u.hourangle, u.deg)).galactic
    l_, b_ = sc.l.deg, sc.b.deg
    IH = float(halpha[hp.ang2pix(NSIDE, l_, b_, lonlat=True)])   # Rayleighs
    EM = max(2.75 * T4 ** 0.9 * IH, 1e-3)                        # pc cm^-6
    DM = c["DM"]
    gold = ("PX" in c and "PXe" in c and c["PX"] > 0 and c["PX"] / c["PXe"] > 3)
    d = 1 / c["PX"] if gold else c.get("DIST_DM", np.nan)
    ne = EM / DM
    s_path = DM * np.log((8000.0) ** 1.5 / max(ne, 1e-6))        # frozen two-phase proxy
    meta.append(dict(DM=DM, d=d, gold=gold, b=b_, EM=EM, ne=ne, s_path=s_path))
R = pd.concat([R, pd.DataFrame(meta)], axis=1).dropna(subset=["DM", "d"])
R["hii"] = (R.EM > R.EM.quantile(0.9)).astype(float)
R.to_csv(os.path.join(DERIVED, "t2_2_sightline_table.csv"), index=False)

# ---------------- statistics ----------------
rng = np.random.default_rng(20260724)
def partial_spearman(y, x, controls):
    def resid_rank(v):
        rk = stats.rankdata(v)
        X = np.column_stack([np.ones(len(v))] + [stats.rankdata(c) for c in controls])
        beta = np.linalg.lstsq(X, rk, rcond=None)[0]
        return rk - X @ beta
    ry, rx = resid_rank(y), resid_rank(x)
    return float(np.corrcoef(ry, rx)[0, 1])

def perm_p(y, x, controls, n=10000):
    obs = partial_spearman(y, x, controls)
    null = np.array([partial_spearman(y, x[rng.permutation(len(x))], controls) for _ in range(n)])
    return obs, float(np.mean(np.abs(null) >= abs(obs)))

PROXIES = ["s_path", "EM", "ne"]                # frozen: primary s_path; secondaries EM, <n_e>
def run_block(df, resp, label):
    ctr = [np.log10(df.DM.values), np.log10(df.d.values)]
    out = {}
    for pr in PROXIES + ["hii"]:
        x = np.log10(np.abs(df[pr].values) + 1e-9) if pr != "hii" else df[pr].values
        r, p = perm_p(df[resp].values, x, ctr)
        out[pr] = dict(partial_rho=round(r, 3), perm_p=round(p, 4))
    pmin = min(v["perm_p"] for k, v in out.items() if k in PROXIES)
    out["_sidak_over_3"] = round(1 - (1 - pmin) ** 3, 4)
    print(f"[{label}] n={len(df)} " + "  ".join(
        f"{k}: rho={v['partial_rho']:+.2f} p={v['perm_p']:.3f}" for k, v in out.items() if k != "_sidak_over_3")
        + f"  | Sidak p={out['_sidak_over_3']:.3f}")
    return out

print("=== T2.2 achromatic response (joint E_RN), all pulsars with covariates ===")
res_all = run_block(R, "logE", "achromatic|all")
res_gold = run_block(R[R.gold], "logE", "achromatic|gold-d")
res_b5 = run_block(R[np.abs(R.b) > 5], "logE", "achromatic|b>5")
print("=== chromatic control (PPTA dm_gp) ===")
Rc = R.dropna(subset=["logE_chrom"])
res_chrom = run_block(Rc, "logE_chrom", "chromatic|PPTA")
# differential (frozen requirement): achromatic must EXCEED chromatic dependence
diff = {pr: round(res_all[pr]["partial_rho"] - res_chrom[pr]["partial_rho"], 3) for pr in PROXIES}
print("differential (achromatic - chromatic) partial rho:", diff)

# loud/quiet logistic (permutation on the primary proxy)
from numpy.linalg import lstsq
def logistic_perm(df, n=5000):
    y = df.loud.values.astype(float)
    X = np.column_stack([np.log10(df.s_path), np.log10(df.DM), np.log10(df.d)])
    X = (X - X.mean(0)) / X.std(0)
    def fit(yv):
        # simple IRLS logistic
        b = np.zeros(X.shape[1] + 1); Xd = np.column_stack([np.ones(len(yv)), X])
        for _ in range(50):
            eta = Xd @ b; mu = 1 / (1 + np.exp(-eta)); W = mu * (1 - mu) + 1e-6
            b = b + lstsq(Xd * W[:, None], (yv - mu), rcond=None)[0]
        return b[1]
    obs = fit(y)
    null = np.array([fit(y[rng.permutation(len(y))]) for _ in range(n)])
    return float(obs), float(np.mean(np.abs(null) >= abs(obs)))
lo, lp = logistic_perm(R)
print(f"logistic loud~s_path|DM,d: beta={lo:+.2f} perm_p={lp:.3f}")

# effect-size gate vs frozen prediction (slope=1 in logE vs log s at fixed DM,d)
m = R[R.gold] if R.gold.sum() >= 10 else R
X = np.column_stack([np.ones(len(m)), np.log10(m.s_path), np.log10(m.DM), np.log10(m.d)])
beta = np.linalg.lstsq(X, m.logE.values, rcond=None)[0]
contrast_pred = float(np.log10(R.s_path.max() / R.s_path.min()))
print(f"measured slope dlogE/dlogs = {beta[1]:+.2f} (frozen prediction: +1.0); "
      f"sample contrast log10(s_max/s_min)={contrast_pred:.2f}")

out = dict(frozen=B4P, table_n=len(R),
           achromatic=dict(all=res_all, gold=res_gold, b_gt5=res_b5),
           chromatic_ppta=res_chrom, differential_rho=diff,
           logistic=dict(beta_s_path=lo, perm_p=lp),
           slope_measured=float(beta[1]), contrast_log10=contrast_pred,
           deviations=["pygedm unavailable (build failed): NE2001 HII tagging replaced by "
                       "EM>90th-pct indicator", "chromatic control PPTA-only (NG uses DMX)"])
json.dump(out, open(os.path.join(DERIVED, "t2_2_results.json"), "w"), indent=1)
print("saved derived/t2_2_results.json + t2_2_sightline_table.csv")

# figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
cmap = {True: "crimson", False: "steelblue"}
ax[0].scatter(R.s_path, 10 ** R.logE, c=R.loud.map(cmap), s=40, edgecolor="k", lw=0.3)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel(r"$s_{\rm path}$ proxy  [pc cm$^{-3}$ ln-units]")
ax[0].set_ylabel(r"achromatic $E_{\rm RN}$ [s$^2$]")
ax[0].set_title(f"T2.2: achromatic excess vs sightline entropy (red=loud)\n"
                f"partial rho={res_all['s_path']['partial_rho']:+.2f}, perm p={res_all['s_path']['perm_p']:.3f}")
ax[0].grid(alpha=0.3, which="both")
labels = PROXIES
xa = [res_all[p]["partial_rho"] for p in labels]
xc = [res_chrom[p]["partial_rho"] for p in labels]
xx = np.arange(len(labels))
ax[1].bar(xx - 0.18, xa, 0.36, label="achromatic (joint E_RN)", color="seagreen")
ax[1].bar(xx + 0.18, xc, 0.36, label="chromatic control (PPTA dm_gp)", color="orange")
ax[1].axhline(0, color="k", lw=1)
ax[1].set_xticks(xx); ax[1].set_xticklabels(labels)
ax[1].set_ylabel(r"partial Spearman $\rho$ (| log DM, log d)")
ax[1].set_title("B4 requires green > orange with green significant — differential test")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "docs/t2_2_entropy_regression.png"), dpi=130)
print("saved docs/t2_2_entropy_regression.png")
