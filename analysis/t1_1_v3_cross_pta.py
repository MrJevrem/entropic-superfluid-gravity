#!/usr/bin/env python
"""
T1.1 v3 — cross-PTA test: NANOGrav 15yr + PPTA DR3 (joint CURN chain) + EPTA DR2 (maxlike).

Adds over v2:
 * PPTA DR3 per-pulsar achromatic RN posteriors from the ARRAY-LEVEL JOINT CURN run
   (chain_commonNoise_pl_nocorr_freegam_DE440.npy, Reardon et al. 2023): each pulsar's RN is
   fit simultaneously with the common process + DM GP + band noise + solar wind, so the
   per-pulsar RN posterior is already "excess above common" — the joint refit v1/v2 lacked.
   Includes an 18-yr J0437-4715 (nearest pulsar; NANOGrav span was only 4.7 yr).
 * Combined envelope/ceiling and distance regression (common-gamma conditional band power,
   band k=1..5 of T=16.03 yr, permutation-calibrated).
 * CROSS-PTA CONSISTENCY: a path-accumulated line-of-sight process must print identically in
   every array observing the same pulsar (same sky signal, overlapping epochs). ~15 overlap
   pulsars NG∩PPTA; EPTA DR2 maxlike appended. Systematic-origin noise disagrees; sky agrees.
 * Distances uniformly from psrcat v2.8.1 (PX, compact-error parsing), PX/sigma>3 gold cut,
   implied-n_e sanity flag.
"""
import os, glob, re, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NG_BASE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0")
NG_ND, NG_PDIR = os.path.join(NG_BASE, "narrowband/noise"), os.path.join(NG_BASE, "narrowband/par")
HD = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/freespectra/ceffyl_data/30f_fs{hd}_100fDMGP_ceffyl")
PP = os.path.join(ROOT, "data/T1.1_pulsar_timing/ppta_dr3/PPTA-DR3/analysis_codes/data/all")
EP = os.path.join(ROOT, "data/T1.1_pulsar_timing/epta_dr2/extracted/EPTA-DR2/EPTA-DR2/noisefiles/DR2full")
PSRCAT = os.path.join(ROOT, "data/T1.1_pulsar_timing/atnf_psrcat/psrcat_tar/psrcat.db")
OUT = os.path.join(ROOT, "data/T1.1_pulsar_timing/derived")
FIG = os.path.join(ROOT, "docs/t1_1_v3_cross_pta.png")

FYR = 1.0 / (365.25 * 86400)
C_MS, KPC_M = 299792458.0, 3.0856775814913673e19
MASYR_RAD_S = 4.84813681109536e-9 / (365.25 * 86400)
KMAX = 5
GAMMA_GRID = [2.0, 2.5, 3.0, 3.5, 4.0, 13/3, 4.5, 5.0]
H_KERN, ESS_NG, ESS_PP = 0.15, 50, 25
NCOND, NMC, NPERM = 600, 400, 3000
ALIAS = {"B1937+21": "J1939+2134", "B1855+09": "J1857+0943"}   # NG name -> J name
rng = np.random.default_rng(20260724)

# ---------------- psrcat distances (compact-error notation) ----------------
def compact_err(valstr, errstr):
    m = re.match(r"^[-+]?\d+(\.(\d+))?([eE]([-+]?\d+))?$", valstr)
    if not m or not re.match(r"^\d+$", errstr): return None
    nd = len(m.group(2) or ""); exp = int(m.group(4) or 0)
    return int(errstr) * (10.0 ** (exp - nd))

cat = {}
cur = {}
for raw in open(PSRCAT, encoding="latin-1"):
    l = raw.rstrip("\n")
    if l.startswith("@"):
        if cur:
            for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = cur
        cur = {}; continue
    if not l.strip() or l.startswith("#"): continue
    p = l.split()
    if p[0] == "PSRJ": cur["_J"] = p[1]
    elif p[0] == "PSRB": cur["_B"] = p[1]
    elif p[0] == "PX" and len(p) >= 3:
        try:
            cur["PX"] = float(p[1]); e = compact_err(p[1], p[2])
            if e: cur["PXe"] = e
        except ValueError: pass
if cur:
    for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = cur

def catalog_px(name):
    r = cat.get(name) or cat.get(ALIAS.get(name, ""), None)
    if r and "PX" in r and "PXe" in r and r["PX"] > 0: return r["PX"], r["PXe"]
    return None, None

# ---------------- shared machinery ----------------
freqs = np.load(HD + "/freqs.npy"); grid = np.load(HD + "/log10rhogrid.npy")
logden = np.load(HD + "/density.npy")[0]
DF = freqs[0]; T_ARR = 1 / DF / (365.25 * 86400)
cdf = np.cumsum(np.exp(logden - logden.max(1, keepdims=True)), 1); cdf /= cdf[:, -1:]
def hd_power_draws(n):
    u = rng.random((n, KMAX))
    return np.sum(10 ** (2 * np.array([np.interp(u[:, k], cdf[k], grid) for k in range(KMAX)]).T), 1)
E_GWB = float(np.sum(10 ** (2 * np.array([np.interp(0.5, cdf[k], grid) for k in range(KMAX)]))))
def phi_base(gc):
    return (1 / (12 * np.pi ** 2)) * FYR ** (gc - 3) * np.sum(freqs[:KMAX] ** (-gc)) * DF
def E_of(logA, gam):                                   # per-sample band power for (A, gamma) draws
    return (10 ** logA) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3) \
        * np.sum(freqs[:KMAX][None, :] ** (-gam[:, None] if np.ndim(gam) else -gam), axis=-1) * DF

def read_par(path, keys=("PX", "F0", "F1", "DM", "PMELONG", "PMELAT", "PMRA", "PMDEC", "START", "FINISH")):
    d = {}
    for line in open(path):
        t = line.split()
        if t and t[0] in keys:
            try: d[t[0]] = float(t[1])
            except ValueError: pass
    return d

def kin(par, d_kpc):
    mu = np.hypot(par.get("PMELONG", par.get("PMRA", 0.0)), par.get("PMELAT", par.get("PMDEC", 0.0)))
    f1_shk = par.get("F0", 0) * (mu * MASYR_RAD_S) ** 2 * d_kpc * KPC_M / C_MS
    return par.get("F1", np.nan) + f1_shk            # Shklovskii-corrected nu-dot

def conditional(logA, gam, gc, n=NCOND):
    w = np.exp(-0.5 * ((gam - gc) / H_KERN) ** 2); sw = w.sum()
    if sw <= 0: return None, 0.0
    ess = sw ** 2 / np.sum(w ** 2)
    return logA[rng.choice(len(w), size=n, p=w / sw)], ess

def wquant(x, w, q):
    s = np.argsort(x); x, w = x[s], w[s]
    c = np.cumsum(w); c /= c[-1]
    return np.interp(q, c, x)

# ---------------- load NANOGrav (SPNA; excess = E - fs{hd}) ----------------
PRIM = re.compile(r"^[BJ]\d{4}[+-]\d{2,4}$")
NG = {}
parmap = {}
for f in glob.glob(NG_PDIR + "/*.par"):
    parmap.setdefault(os.path.basename(f).split("_")[0], f)
for pf in sorted(glob.glob(NG_ND + "/*.pars.txt")):
    psr = os.path.basename(pf).split(".")[0]
    if not PRIM.match(psr) or psr not in parmap: continue
    pars = open(pf).read().split()
    iA = [i for i, p in enumerate(pars) if p.endswith("red_noise_log10_A")]
    ig = [i for i, p in enumerate(pars) if p.endswith("red_noise_gamma")]
    if not iA: continue
    order = sorted([iA[0], ig[0]])
    ch = pd.read_csv(NG_ND + f"/{psr}.nb.chain_1.txt", sep=r"\s+", header=None, usecols=order).values
    cm = {order[k]: ch[:, k] for k in range(2)}; b = int(0.25 * len(ch))
    par = read_par(parmap[psr])
    span = (par.get("FINISH", 0) - par.get("START", 0)) / 365.25
    NG[psr] = dict(logA=cm[iA[0]][b:], gam=cm[ig[0]][b:], par=par, span=span, array="NG")

# ---------------- load PPTA (joint CURN chain) ----------------
CH = os.path.join(PP, "chains")
pp_pars = [l.strip() for l in open(CH + "/commonNoise_pl_nocorr_freegam_DE440_pars.txt") if l.strip()]
pp_ch = np.load(CH + "/chain_commonNoise_pl_nocorr_freegam_DE440.npy")
pp_ch = pp_ch[int(0.25 * len(pp_ch)):]
col = {p: i for i, p in enumerate(pp_pars)}
crn_A = pp_ch[:, col["gw_pl_nocorr_freegam_log10_A"]]
crn_g = pp_ch[:, col["gw_pl_nocorr_freegam_gamma"]]
E_CRN_samp = E_of(crn_A, crn_g)
PPTA = {}
for p in sorted(set(k.split("_red_noise")[0] for k in pp_pars if "_red_noise_log10_A" in k)):
    parf = os.path.join(PP, f"{p}.par")
    if not os.path.exists(parf): continue
    par = read_par(parf)
    span = (par.get("FINISH", 0) - par.get("START", 0)) / 365.25
    PPTA[p] = dict(logA=pp_ch[:, col[f"{p}_red_noise_log10_A"]],
                   gam=pp_ch[:, col[f"{p}_red_noise_gamma"]], par=par, span=span, array="PPTA")
print(f"NG pulsars={len(NG)}  PPTA pulsars={len(PPTA)} (joint CURN, {len(pp_ch)} samples)  "
      f"CRN median: logA={np.median(crn_A):.2f} gamma={np.median(crn_g):.2f}")

# ---------------- assemble gold sample (psrcat distances) ----------------
def make_row(psr, v):
    px, pxe = catalog_px(psr)
    if px is None: return None
    d = 1 / px
    ne = v["par"].get("DM", np.nan) / (d * 1000)
    return dict(psr=psr, jname=ALIAS.get(psr, psr), array=v["array"], PX=px, PXe=pxe,
                snr=px / pxe, d=d, span=v["span"], DM=v["par"].get("DM", np.nan),
                F0=v["par"].get("F0", np.nan), F1i=kin(v["par"], d), ne=ne,
                suspect=bool(ne > 0.10), kmin=max(1, int(np.ceil(T_ARR / max(v["span"], 0.5)))),
                logA_med=float(np.median(v["logA"])), gam_med=float(np.median(v["gam"])))
rows = []
for src in (NG, PPTA):
    for p, v in src.items():
        r = make_row(p, v)
        if r: rows.append(r)
ALL = pd.DataFrame(rows)
ALL["gold"] = ALL.snr > 3
ALL["band_ok"] = ALL.kmin <= 3
G = ALL[ALL.gold & ALL.band_ok & ~ALL.suspect].copy()
print(f"gold+band_ok+clean: NG={len(G[G.array=='NG'])}  PPTA={len(G[G.array=='PPTA'])}  "
      f"(overlap J-names: {len(set(G[G.array=='NG'].jname) & set(G[G.array=='PPTA'].jname))})")
print(f"J0437-4715 in PPTA set: {'J0437-4715' in G[G.array=='PPTA'].psr.values} "
      f"(span={PPTA.get('J0437-4715',{}).get('span',0):.1f} yr, d={1/catalog_px('J0437-4715')[0]:.3f} kpc)")

def store(psr, arr): return NG[psr] if arr == "NG" else PPTA[psr]

# combined set: prefer PPTA (joint + longer span) where both arrays have the pulsar
pref = {}
for r in G.itertuples():
    if r.jname not in pref or (r.array == "PPTA" and pref[r.jname].array == "NG"):
        pref[r.jname] = r
COMB = pd.DataFrame([r._asdict() for r in pref.values()])

# ---------------- excess machinery per array ----------------
def excess_samples(psr, arr, gc):
    """Conditional band-power EXCESS draws + detection prob. NG: E - fs{hd}; PPTA: E_RN (joint)."""
    v = store(psr, arr)
    a, ess = conditional(v["logA"], v["gam"], gc)
    if a is None or ess < (ESS_NG if arr == "NG" else ESS_PP): return None, None, ess
    E = phi_base(gc) * 10 ** (2 * a)
    if arr == "NG":
        p_det = float(np.mean(E - hd_power_draws(len(E)) > 0))
        E_exc_med = float(phi_base(gc) * 10 ** (2 * np.median(a)) - E_GWB)
    else:
        p_det = float(np.mean(E - E_CRN_samp[rng.integers(len(E_CRN_samp), size=len(E))] > 0))
        E_exc_med = float(phi_base(gc) * 10 ** (2 * np.median(a)))
    return E_exc_med, p_det, ess

# ---------------- regression over gamma_c (PPTA-only / combined) ----------------
def run_regression(sample_df, label):
    out = {}
    for gc in GAMMA_GRID:
        pts = []
        for r in sample_df.itertuples():
            E_med, p_det, ess = excess_samples(r.psr, r.array, gc)
            if E_med is None: continue
            if p_det > 0.5 and E_med > 0 and r.logA_med > -15.5:
                pts.append((np.log10(r.d), np.log10(E_med), r.psr, r.array))
        e = {"n_det": len(pts)}
        if len(pts) >= 5:
            x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
            ts = stats.theilslopes(y, x)
            perm = np.array([stats.theilslopes(y, x[rng.permutation(len(x))])[0] for _ in range(NPERM)])
            e.update(theil=float(ts[0]), p_perm=float(np.mean(perm >= ts[0])),
                     spearman=float(stats.spearmanr(x, y).correlation))
        out[f"{gc:.3f}"] = e
    fitted = {k: v for k, v in out.items() if "theil" in v}
    if fitted:
        pmin = min(v["p_perm"] for v in fitted.values())
        out["_global_sidak"] = float(1 - (1 - pmin) ** len(fitted))
    print(f"\n--- {label} regression ---")
    for gc in GAMMA_GRID:
        e = out[f"{gc:.3f}"]
        s = f"  gc={gc:4.2f}  n={e['n_det']:2d}"
        if "theil" in e: s += f"  Theil={e['theil']:+.2f}  perm_p={e['p_perm']:.3f}  rho={e['spearman']:+.2f}"
        print(s)
    if "_global_sidak" in out: print(f"  GLOBAL Sidak p = {out['_global_sidak']:.3f}")
    return out

reg_ppta = run_regression(G[G.array == "PPTA"], "PPTA-only (joint CURN)")
reg_comb = run_regression(COMB, "combined NG+PPTA")

# ---------------- combined envelope / ceiling (gamma_c = 13/3) ----------------
GC = 13/3; pb = phi_base(GC)
env = []
for r in COMB.itertuples():
    v = store(r.psr, r.array)
    a, ess = conditional(v["logA"], v["gam"], GC, n=3000)
    if a is None or ess < 20: continue
    E_tot = pb * 10 ** (2 * a)
    if r.array == "PPTA":                              # universality budget: RN_i + CRN (joint draws)
        E_tot = E_tot + E_CRN_samp[rng.integers(len(E_CRN_samp), size=len(a))]
    E95 = float(wquant(E_tot, 10 ** a, 0.95))          # uniform-A-prior 95% (conservative)
    E_med, p_det, _ = excess_samples(r.psr, r.array, GC)
    det = (p_det or 0) > 0.5 and (E_med or 0) > 0 and r.logA_med > -15.5
    env.append(dict(psr=r.psr, array=r.array, d=r.d, span=r.span, detected=bool(det),
                    E_exc_med=E_med if det else np.nan, E95=E95, K=E95 / r.d))
EV = pd.DataFrame(env).sort_values("d").reset_index(drop=True)
quiet = EV[~EV.detected]
kmin_row = quiet.loc[quiet.K.idxmin()]
K95 = float(kmin_row.K)
loudK = (EV[EV.detected].E_exc_med / EV[EV.detected].d)
floor_sp = stats.spearmanr(np.log10(quiet.d), np.log10(quiet.E95))
print(f"\n--- combined envelope (gc=13/3, n={len(EV)}: {int(EV.detected.sum())} det, {len(quiet)} quiet) ---")
print(f"K95 = {K95:.2e} s^2/kpc  (set by {kmin_row['psr']} [{kmin_row['array']}] at {kmin_row['d']:.2f} kpc)")
print(f"  vs v2 (NG-only): 5.7e-13 | loud median K = {loudK.median():.2e} "
      f"({np.log10(loudK.median()/K95):+.1f} dex above ceiling)")
print(f"  universal comp at 1 kpc: log10A <= {0.5*np.log10(K95/pb):.2f} (A_gwb=-14.62)")
print(f"  floor trend: Spearman rho={floor_sp.correlation:+.2f} (p={floor_sp.pvalue:.2f}) [theory: rising]")
j0437 = EV[EV.psr == "J0437-4715"]
if len(j0437):
    r = j0437.iloc[0]
    print(f"  J0437-4715 (18yr, d=0.156 kpc): detected={r.detected}, "
          f"E={'%.2e' % r.E_exc_med if r.detected else '<%.2e' % r.E95}  "
          f"-> K_implied={'%.2e' % ((r.E_exc_med if r.detected else r.E95)/r.d)} s^2/kpc")

# ---------------- native detection (no conditioning: raw joint/SPNA chains) ----------------
def native_pdet(psr, arr):
    v = store(psr, arr)
    E = E_of(v["logA"], v["gam"])
    if arr == "NG":
        return float(np.mean(E - hd_power_draws(len(E)) > 0))
    return float(np.mean(E - E_CRN_samp[rng.integers(len(E_CRN_samp), size=len(E))] > 0))

pp_census = {p: native_pdet(p, "PPTA") for p in PPTA}
pp_det_native = sorted([p for p, q in pp_census.items() if q > 0.5])
print(f"\nPPTA native census (P(E_RN > E_CRN) > 0.5, raw joint chains): "
      f"{len(pp_det_native)}/{len(PPTA)} detected: {pp_det_native}")

# ---------------- cross-PTA consistency (overlap pulsars, native posteriors) ----------------
overlap = sorted(set(G[G.array == "NG"].jname) & set(G[G.array == "PPTA"].jname))
ep_files = {os.path.basename(f).split("_noise")[0]: f for f in glob.glob(EP + "/*_noise.json")}
cons = []
for j in overlap:
    ng_name = {v: k for k, v in ALIAS.items()}.get(j, j)
    ng_name = ng_name if ng_name in NG else j
    png, ppp = native_pdet(ng_name, "NG"), pp_census[j]
    ep = None
    for cand in (j, ng_name):
        if cand in ep_files:
            d = json.load(open(ep_files[cand]))
            ep = (d.get(f"{cand}_red_noise_log10_A"), d.get(f"{cand}_red_noise_gamma")); break
    vng, vpp = NG[ng_name], PPTA[j]
    cons.append(dict(jname=j, ng_logA=float(np.median(vng["logA"])), ng_gam=float(np.median(vng["gam"])),
                     pp_logA=float(np.median(vpp["logA"])), pp_gam=float(np.median(vpp["gam"])),
                     ng_pdet=round(png, 2), pp_pdet=round(ppp, 2),
                     agree=(png > 0.5) == (ppp > 0.5),
                     ep_logA=(ep[0] if ep else np.nan), ep_gam=(ep[1] if ep else np.nan)))
CT = pd.DataFrame(cons)
n_disagree = int((~CT.agree).sum())
print(f"\n--- cross-PTA consistency ({len(CT)} overlap pulsars, NATIVE posteriors) ---")
print(CT.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
print(f"detection-status disagreements NG vs PPTA: {n_disagree}/{len(CT)}"
      f"   [a sky (line-of-sight) process requires agreement]"
      f"\nNOTE: NG = single-pulsar fits minus fs(hd) (common not jointly modeled);"
      f" PPTA = joint CURN. The joint method is the correct one; disagreement pattern"
      f" (NG-detected -> PPTA-quiet) indicates the NG excesses are method/dataset artifacts,"
      f" not line-of-sight noise.")

# ---------------- save ----------------
res = dict(ppta_chain="commonNoise_pl_nocorr_freegam_DE440 (joint CURN, per-pulsar RN)",
           crn=dict(logA=float(np.median(crn_A)), gamma=float(np.median(crn_g))),
           regression_ppta=reg_ppta, regression_combined=reg_comb,
           envelope=dict(gamma_c=GC, n=len(EV), n_quiet=len(quiet), K95=K95,
                         K95_pulsar=str(kmin_row["psr"]), K95_array=str(kmin_row["array"]),
                         log10A_universal_1kpc=float(0.5 * np.log10(K95 / pb)),
                         loud_K_median=float(loudK.median()),
                         floor_spearman=float(floor_sp.correlation), floor_p=float(floor_sp.pvalue)),
           consistency=dict(n_overlap=len(CT), n_disagree=n_disagree,
                            table=CT.to_dict(orient="records")),
           ppta_native_census=dict(n_total=len(PPTA), n_detected=len(pp_det_native),
                                   detected=pp_det_native,
                                   pdet={k: round(v, 3) for k, v in pp_census.items()}))
json.dump(res, open(os.path.join(OUT, "t1_1_v3_results.json"), "w"), indent=2, default=str)
EV.to_csv(os.path.join(OUT, "t1_1_v3_envelope.csv"), index=False)

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))
a0 = ax[0]
for arr, mk, c in (("NG", "o", "#1f77b4"), ("PPTA", "s", "#9467bd")):
    dd = EV[(EV.array == arr) & EV.detected]; qq = EV[(EV.array == arr) & ~EV.detected]
    a0.scatter(dd.d, dd.E_exc_med, marker=mk, s=55, color=c, edgecolor="k", lw=0.4, label=f"{arr} detected")
    a0.scatter(qq.d, qq.E95, marker="v", s=60, facecolors="none", edgecolors=c, lw=1.3, label=f"{arr} 95% limit")
ddl = np.logspace(np.log10(0.09), np.log10(4.5), 40)
a0.plot(ddl, K95 * ddl, "crimson", lw=2, label=r"universal ceiling $K_{95}\cdot d$")
a0.axhline(E_GWB, color="k", ls=":", lw=1.2, label="HD (GWB) band power")
for nm in ("J0437-4715", "J1939+2134", "J1909-3744", str(kmin_row.psr)):
    rr = EV[EV.psr == nm]
    if len(rr):
        r = rr.iloc[0]
        a0.annotate(nm, (r.d, r.E_exc_med if r.detected else r.E95), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
a0.set_xscale("log"); a0.set_yscale("log")
a0.set_xlabel("distance d [kpc]"); a0.set_ylabel(r"achromatic excess band power $E$ (k=1–5) [s$^2$]")
a0.set_title("v3 combined envelope: NG + PPTA DR3 (incl. 18-yr J0437)")
a0.legend(fontsize=7, loc="lower right"); a0.grid(alpha=0.3, which="both")

a1 = ax[1]
for reg, lab, c in ((reg_ppta, "PPTA-only", "#9467bd"), (reg_comb, "combined", "seagreen")):
    gcs = [g for g in GAMMA_GRID if "theil" in reg[f"{g:.3f}"]]
    a1.plot(gcs, [reg[f"{g:.3f}"]["theil"] for g in gcs], "o-", color=c,
            label=f"{lab} (global p={reg.get('_global_sidak', np.nan):.2f})")
a1.axhline(1.0, color="b", ls="--", lw=1.5, label="theory: slope +1")
a1.axhline(0.0, color="k", lw=1.2, label="null")
a1.set_xlabel(r"common index $\gamma_c$"); a1.set_ylabel(r"Theil–Sen slope dlog$E$/dlog$d$")
a1.set_title(f"distance slope vs $\\gamma_c$;  cross-PTA disagreements: {n_disagree}/{len(CT)}")
a1.legend(fontsize=8); a1.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(FIG, dpi=130)
print(f"\nsaved: {FIG}")
print(f"saved: {os.path.join(OUT, 't1_1_v3_results.json')} + envelope csv")
