#!/usr/bin/env python
"""
T1.1 v4 — joint-vs-joint: NANOGrav 15yr CURN chain (la_forge core) vs PPTA DR3 CURN chain.

Both sides now use per-pulsar achromatic RN fit JOINTLY with the common process, inside each
array's likelihood ("excess above common" is symmetric and method-matched):
  NG:   curn_14f_pl_vg.core  (67 psr x [gamma, log10_A] + gw_gamma/gw_log10_A; 42k samples)
  PPTA: chain_commonNoise_pl_nocorr_freegam_DE440.npy (30 psr + CRN)

Products:
 1. NG-joint native census: P(E_RN,i > E_gw) rowwise (paired samples) — parallel to PPTA's 2/30.
 2. SPNA -> joint "dissolution": which v1/v3 NG excesses survive joint modeling.
 3. Joint-vs-joint consistency table for overlap pulsars (+ EPTA maxlike column).
 4. Fully-joint combined envelope ceiling K95 (budget = E_RN + E_common each side).
"""
import os, glob, re, json
import numpy as np
import pandas as pd
import h5py
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NG_BASE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0")
NG_ND, NG_PDIR = os.path.join(NG_BASE, "narrowband/noise"), os.path.join(NG_BASE, "narrowband/par")
CORE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/gwb_cores/curn_14f_pl_vg.core")
HD_CORE = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/gwb_cores/hd_14f_pl_vg.core")
FS = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/freespectra/ceffyl_data/30f_fs{hd}_100fDMGP_ceffyl")
PP = os.path.join(ROOT, "data/T1.1_pulsar_timing/ppta_dr3/PPTA-DR3/analysis_codes/data/all")
EP = os.path.join(ROOT, "data/T1.1_pulsar_timing/epta_dr2/extracted/EPTA-DR2/EPTA-DR2/noisefiles/DR2full")
PSRCAT = os.path.join(ROOT, "data/T1.1_pulsar_timing/atnf_psrcat/psrcat_tar/psrcat.db")
OUT = os.path.join(ROOT, "data/T1.1_pulsar_timing/derived")

FYR = 1.0 / (365.25 * 86400)
KMAX = 5
ALIAS = {"B1937+21": "J1939+2134", "B1855+09": "J1857+0943"}
rng = np.random.default_rng(20260724)

freqs = np.load(FS + "/freqs.npy")[:KMAX]
DF = np.load(FS + "/freqs.npy")[0]
T_ARR = 1 / DF / (365.25 * 86400)
def E_of(logA, gam):
    gam = np.atleast_1d(gam)
    return ((10 ** np.atleast_1d(logA)) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3)
            * np.sum(freqs[None, :] ** (-gam[:, None]), axis=1) * DF)
def phi_base(gc):
    return (1 / (12 * np.pi ** 2)) * FYR ** (gc - 3) * np.sum(freqs ** (-gc)) * DF

# ---------------- NG joint core ----------------
with h5py.File(CORE) as f:
    pars = [p.decode() if isinstance(p, bytes) else str(p) for p in f["params"][:]]
    burn = int(f["metadata/burn"][()])
    ch = f["chain"][:]
ch = ch[burn:]
col = {p: i for i, p in enumerate(pars)}
E_gw = E_of(ch[:, col["gw_log10_A"]], ch[:, col["gw_gamma"]])
NGJ = {}
for p in sorted(set(k.split("_red_noise")[0] for k in pars if "_red_noise_log10_A" in k)):
    NGJ[p] = dict(logA=ch[:, col[f"{p}_red_noise_log10_A"]], gam=ch[:, col[f"{p}_red_noise_gamma"]])
gwA, gwG = np.median(ch[:, col["gw_log10_A"]]), np.median(ch[:, col["gw_gamma"]])
print(f"NG joint CURN core: {len(NGJ)} pulsars, {len(ch)} post-burn samples; "
      f"gw: logA={gwA:.2f} gamma={gwG:.2f}")

with h5py.File(HD_CORE) as f:
    hpars = [p.decode() if isinstance(p, bytes) else str(p) for p in f["params"][:]]
    hburn = int(f["metadata/burn"][()])
    hch = f["chain"][:][hburn:]
hcol = {p: i for i, p in enumerate(hpars)}
hd_A = next(p for p in hpars if re.match(r"^gw.*log10_A$", p))
hd_G = next(p for p in hpars if re.match(r"^gw.*gamma$", p))
E_gw_hd = E_of(hch[:, hcol[hd_A]], hch[:, hcol[hd_G]])

def pdet_joint(store, Ecom, p):
    E = E_of(store[p]["logA"], store[p]["gam"])
    return float(np.mean(E > Ecom))            # rowwise-paired: same chain rows

# ---------------- NG census (CURN + HD-common robustness) ----------------
census = {p: pdet_joint(NGJ, E_gw, p) for p in NGJ}
census_hd = {p: float(np.mean(E_of(hch[:, hcol[f"{p}_red_noise_log10_A"]],
                                   hch[:, hcol[f"{p}_red_noise_gamma"]]) > E_gw_hd))
             for p in NGJ if f"{p}_red_noise_log10_A" in hcol}
det = sorted([p for p, q in census.items() if q > 0.5])
det_hd = sorted([p for p, q in census_hd.items() if q > 0.5])
print(f"\nNG-joint census (P(E_RN > E_gw) > 0.5): {len(det)}/{len(NGJ)} detected")
print("  CURN-common:", det)
print(f"  HD-common robustness: {len(det_hd)}/{len(census_hd)}:", det_hd)

# ---------------- SPNA -> joint dissolution ----------------
PRIM = re.compile(r"^[BJ]\d{4}[+-]\d{2,4}$")
grid = np.load(FS + "/log10rhogrid.npy"); logden = np.load(FS + "/density.npy")[0]
cdf = np.cumsum(np.exp(logden - logden.max(1, keepdims=True)), 1); cdf /= cdf[:, -1:]
def hd_power_draws(n):
    u = rng.random((n, KMAX))
    return np.sum(10 ** (2 * np.array([np.interp(u[:, k], cdf[k], grid) for k in range(KMAX)]).T), 1)
spna = {}
for pf in sorted(glob.glob(NG_ND + "/*.pars.txt")):
    psr = os.path.basename(pf).split(".")[0]
    if not PRIM.match(psr) or psr not in NGJ: continue
    pr = open(pf).read().split()
    iA = [i for i, p in enumerate(pr) if p.endswith("red_noise_log10_A")]
    ig = [i for i, p in enumerate(pr) if p.endswith("red_noise_gamma")]
    if not iA: continue
    order = sorted([iA[0], ig[0]])
    c = pd.read_csv(NG_ND + f"/{psr}.nb.chain_1.txt", sep=r"\s+", header=None, usecols=order).values
    cm = {order[k]: c[:, k] for k in range(2)}; b = int(0.25 * len(c))
    E = E_of(cm[iA[0]][b:], cm[ig[0]][b:])
    spna[psr] = float(np.mean(E - hd_power_draws(len(E)) > 0))
spna_det = {p for p, q in spna.items() if q > 0.5}
dissolved = sorted(spna_det - set(det))
survived = sorted(spna_det & set(det))
print(f"\nSPNA-detected (v3 method): {len(spna_det)} | survive joint modeling: {len(survived)} "
      f"| dissolve: {len(dissolved)}")
print("  survived: ", survived)
print("  dissolved:", dissolved)

# ---------------- PPTA joint (as v3) ----------------
CH = os.path.join(PP, "chains")
pp_pars = [l.strip() for l in open(CH + "/commonNoise_pl_nocorr_freegam_DE440_pars.txt") if l.strip()]
pp_ch = np.load(CH + "/chain_commonNoise_pl_nocorr_freegam_DE440.npy")
pp_ch = pp_ch[int(0.25 * len(pp_ch)):]
pcol = {p: i for i, p in enumerate(pp_pars)}
E_crn = E_of(pp_ch[:, pcol["gw_pl_nocorr_freegam_log10_A"]], pp_ch[:, pcol["gw_pl_nocorr_freegam_gamma"]])
PPJ = {p: dict(logA=pp_ch[:, pcol[f"{p}_red_noise_log10_A"]], gam=pp_ch[:, pcol[f"{p}_red_noise_gamma"]])
       for p in sorted(set(k.split("_red_noise")[0] for k in pp_pars if "_red_noise_log10_A" in k))}
pp_census = {p: pdet_joint(PPJ, E_crn, p) for p in PPJ}

# ---------------- joint-vs-joint consistency table ----------------
ep_files = {os.path.basename(f).split("_noise")[0]: f for f in glob.glob(EP + "/*_noise.json")}
J2NG = {v: k for k, v in ALIAS.items()}
overlap = sorted(set(ALIAS.get(p, p) for p in NGJ) & set(PPJ))
rows = []
for j in overlap:
    ng = J2NG.get(j, j)
    ep = None
    for cand in (j, ng):
        if cand in ep_files:
            d = json.load(open(ep_files[cand]))
            ep = d.get(f"{cand}_red_noise_log10_A"); break
    q_ng, q_pp = census[ng], pp_census[j]
    rows.append(dict(jname=j,
                     ngJ_logA=float(np.median(NGJ[ng]["logA"])), ngJ_gam=float(np.median(NGJ[ng]["gam"])),
                     ngJ_pdet=round(q_ng, 2),
                     pp_logA=float(np.median(PPJ[j]["logA"])), pp_gam=float(np.median(PPJ[j]["gam"])),
                     pp_pdet=round(q_pp, 2),
                     agree=(q_ng > 0.5) == (q_pp > 0.5),
                     spna_pdet=round(spna.get(ng, np.nan), 2), ep_logA=ep if ep is not None else np.nan))
CT = pd.DataFrame(rows)
n_dis = int((~CT.agree).sum())
print(f"\n--- JOINT-vs-JOINT consistency ({len(CT)} overlap pulsars) ---")
print(CT.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
print(f"disagreements: {n_dis}/{len(CT)}   (v3 SPNA-vs-joint had 6/12)")

# ---------------- fully-joint combined envelope (gamma_c = 13/3) ----------------
def compact_err(valstr, errstr):
    m = re.match(r"^[-+]?\d+(\.(\d+))?([eE]([-+]?\d+))?$", valstr)
    if not m or not re.match(r"^\d+$", errstr): return None
    nd = len(m.group(2) or ""); exp = int(m.group(4) or 0)
    return int(errstr) * (10.0 ** (exp - nd))
cat, cur = {}, {}
for raw in open(PSRCAT, encoding="latin-1"):
    l = raw.rstrip("\n")
    if l.startswith("@"):
        if cur:
            for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = cur
        cur = {}; continue
    if not l.strip() or l.startswith("#"): continue
    t = l.split()
    if t[0] == "PSRJ": cur["_J"] = t[1]
    elif t[0] == "PSRB": cur["_B"] = t[1]
    elif t[0] == "PX" and len(t) >= 3:
        try:
            cur["PX"] = float(t[1]); e = compact_err(t[1], t[2])
            if e: cur["PXe"] = e
        except ValueError: pass
if cur:
    for nm in filter(None, [cur.get("_J"), cur.get("_B")]): cat[nm] = cur

def px_of(name):
    r = cat.get(name) or cat.get(ALIAS.get(name, ""), None)
    return (r["PX"], r["PXe"]) if r and "PX" in r and "PXe" in r and r["PX"] > 0 else (None, None)

def par_meta(path):
    d = {}
    for line in open(path):
        t = line.split()
        if t and t[0] in ("DM", "START", "FINISH"):
            try: d[t[0]] = float(t[1])
            except ValueError: pass
    return d

def conditional(logA, gam, gc, n=3000, h=0.15):
    w = np.exp(-0.5 * ((gam - gc) / h) ** 2); sw = w.sum()
    if sw <= 0: return None
    return logA[rng.choice(len(w), size=n, p=w / sw)]
def wquant(x, w, q):
    s = np.argsort(x); x, w = x[s], w[s]
    c = np.cumsum(w); c /= c[-1]
    return np.interp(q, c, x)

GC = 13 / 3; pb = phi_base(GC)
parmap = {}
for f in glob.glob(NG_PDIR + "/*.par"):
    parmap.setdefault(os.path.basename(f).split("_")[0], f)
env = []
for arr, store, Ecom, cenmap in (("NG", NGJ, E_gw, census), ("PPTA", PPJ, E_crn, pp_census)):
    for p in store:
        px, pxe = px_of(p)
        if px is None or px / pxe <= 3: continue
        pf = parmap.get(p) if arr == "NG" else os.path.join(PP, f"{p}.par")
        if not pf or not os.path.exists(pf): continue
        m = par_meta(pf)
        span = (m.get("FINISH", 0) - m.get("START", 0)) / 365.25
        if max(1, int(np.ceil(T_ARR / max(span, 0.5)))) > 3: continue          # band validity
        d = 1 / px
        if (m.get("DM", 0) / (d * 1000)) > 0.10: continue                       # suspect PX
        a = conditional(store[p]["logA"], store[p]["gam"], GC)
        if a is None: continue
        E_tot = pb * 10 ** (2 * a) + Ecom[rng.integers(len(Ecom), size=len(a))]
        env.append(dict(psr=p, arr=arr, jname=ALIAS.get(p, p), d=d,
                        detected=cenmap[p] > 0.5,
                        E95=float(wquant(E_tot, 10 ** a, 0.95))))
EV = pd.DataFrame(env)
EV = EV.loc[EV.groupby("jname").arr.idxmax()] if False else EV                  # keep both arrays
quiet = EV[~EV.detected].copy(); quiet["K"] = quiet.E95 / quiet.d
kr = quiet.loc[quiet.K.idxmin()]
K95 = float(kr.K)
print(f"\n--- fully-joint envelope (gc=13/3, {len(EV)} entries, {len(quiet)} quiet) ---")
print(f"K95 = {K95:.2e} s^2/kpc  (set by {kr['psr']} [{kr['arr']}] at {kr['d']:.2f} kpc)")
print(f"  universal comp at 1 kpc: log10A <= {0.5*np.log10(K95/pb):.2f}")
print(f"  (v2 NG-SPNA: 5.72e-13 J1640+2224 | v3 mixed: 5.74e-13 J1909-3744 [PPTA])")

json.dump(dict(ng_core="curn_14f_pl_vg.core (67 psr joint CURN, 14f common)",
               gw=dict(logA=float(gwA), gamma=float(gwG)),
               ng_census=dict(n=len(NGJ), n_det=len(det), detected=det,
                              hd_common_n_det=len(det_hd), hd_common_detected=det_hd,
                              pdet={k: round(v, 3) for k, v in census.items()}),
               dissolution=dict(spna_detected=sorted(spna_det), survived=survived, dissolved=dissolved),
               consistency=dict(n=len(CT), n_disagree=n_dis, table=CT.to_dict(orient="records")),
               envelope=dict(K95=K95, pulsar=str(kr["psr"]), array=str(kr["arr"]),
                             log10A_1kpc=float(0.5 * np.log10(K95 / pb)))),
          open(os.path.join(OUT, "t1_1_v4_results.json"), "w"), indent=2, default=str)
CT.to_csv(os.path.join(OUT, "t1_1_v4_consistency.csv"), index=False)
print(f"\nsaved: {os.path.join(OUT,'t1_1_v4_results.json')} + consistency csv")
