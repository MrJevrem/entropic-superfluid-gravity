#!/usr/bin/env python
"""
T3.7 — the mixed-method sample (Paper V sec.7 upgrade (i); OPEN_ITEMS E-4 slice).

Design: g_obs from CAUSTIC masses (HeCS-SZ, galaxy-dynamics — independent of the
X-ray hydrostatic chain), sigma from member dispersions (same catalog), g_bar from
ACCEPT entropy profiles + <kT> (n_e = (kT/K)^{3/2}) with the frozen nominal-BCG
convention. This breaks the mechanical entanglement flagged in Paper V sec.4 and
enables the deferred mass-at-fixed-sigma test. Forks: SZ-mass substitute for the
caustic mass; Duffy c(M,z) concentration (baseline) vs c=5 (flat).
Statistic and window are the FROZEN ones (T3.0): D = <log10 g_obs - log10 g_RAR>
over g_bar in [1e-12, 1e-10]; a0 = 1.0801e-10; >=3 window points.
"""
import os, sys, json, re, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import enforce, DERIVED, ROOT

G = 6.674e-11; MSUN = 1.989e30; KPC = 3.0857e19; A0 = 1.0801e-10
H0 = 70/3.086e19; OM, OL = 0.3, 0.7
gRAR = lambda gb: gb/(1 - np.exp(-np.sqrt(gb/A0)))

MDIR = os.path.join(ROOT, "data/T3/mixed"); os.makedirs(MDIR, exist_ok=True)
HECS = os.path.join(MDIR, "hecs_sz_table4.dat")
if not os.path.exists(HECS):
    urllib.request.urlretrieve("https://cdsarc.cds.unistra.fr/ftp/J/ApJ/819/63/table4.dat", HECS)
    print("staged HeCS-SZ table4")
enforce(staged=[HECS], allow_dirty_paths=("t3_7_mixed_method.py",))

# ---- parse HeCS-SZ table4 (caustic masses + sigma + SZ masses) ----
CL = []
for L in open(HECS):
    try:
        cid = L[0:12].strip(); z = float(L[33:39]); sig = float(L[40:44])
        m200 = float(L[52:57])*1e14; msz = L[63:68].strip()
        CL.append(dict(id=cid, z=z, sig=sig, M=m200, MSZ=float(msz)*1e14 if msz else None))
    except Exception: continue
print(f"HeCS-SZ clusters parsed: {len(CL)} (sigma {min(c['sig'] for c in CL):.0f}-{max(c['sig'] for c in CL):.0f})")

# ---- parse ACCEPT (staged pre-existing): table1 <kT>, z; table5 K0,K100,alpha ----
ACC = {}
for L in open(os.path.join(ROOT, "data/T3/accept/table1.dat")):
    try:
        nm = L[0:18].strip(); z = float(L[60:66]); kt = float(L[67:72])
        if nm not in ACC: ACC[nm] = dict(z=z, kT=kt)
    except Exception: continue
for L in open(os.path.join(ROOT, "data/T3/accept/table5.dat")):
    try:
        nm = L[0:18].strip()
        if nm in ACC and "K0" not in ACC[nm]:
            ACC[nm].update(K0=float(L[32:37]), K100=float(L[51:58]), alpha=float(L[66:70]))
    except Exception: continue
ACC = {k: v for k, v in ACC.items() if "K0" in v}
print(f"ACCEPT clusters with entropy fits: {len(ACC)}")

norm = lambda s: (lambda m: f"A{int(m.group(1))}" if m else re.sub(r"\s+", "", s.upper()))(
    re.search(r"(?:ABELL|^A)\s*0*(\d+)", s.upper()))
ACCN = {norm(k): v for k, v in ACC.items()}

def rho_c(z): return 3*(H0**2*(OM*(1+z)**3 + OL))/(8*np.pi*G)
def nfw_g(r, M200, z, cfork="duffy"):
    r200 = (3*M200*MSUN/(4*np.pi*200*rho_c(z)))**(1/3)
    c = 5.0 if cfork == "flat" else 5.71*(M200/(2e12/0.7))**-0.084*(1+z)**-0.47
    f = lambda x: np.log(1+x) - x/(1+x)
    return G*M200*MSUN*f(r*c/r200)/f(c)/r**2

rows, skipped = [], 0
for c in CL:
    a = ACCN.get(norm(c["id"]))
    if not a: continue
    r = np.linspace(20, 2000, 400)*KPC
    K = a["K0"] + a["K100"]*(r/(100*KPC))**a["alpha"]
    ne = np.maximum(a["kT"]/np.maximum(K, 1e-3), 0)**1.5          # cm^-3
    rho_g = 1.97e-21*ne                                            # kg/m^3
    Mg = np.concatenate([[0], np.cumsum(4*np.pi*r[1:]**2*rho_g[1:]*np.diff(r))])
    gbar = G*(Mg + 1e12*MSUN)/r**2
    w = (gbar >= 1e-12) & (gbar <= 1e-10)
    if w.sum() < 3: skipped += 1; continue
    D = float(np.mean(np.log10(nfw_g(r[w], c["M"], c["z"])) - np.log10(gRAR(gbar[w]))))
    Dc5 = float(np.mean(np.log10(nfw_g(r[w], c["M"], c["z"], "flat")) - np.log10(gRAR(gbar[w]))))
    Dsz = None
    if c["MSZ"]:
        Dsz = float(np.mean(np.log10(nfw_g(r[w], c["MSZ"], c["z"])) - np.log10(gRAR(gbar[w]))))
    rows.append(dict(id=c["id"], z=c["z"], sig=c["sig"], M=c["M"], D=D, D_c5=Dc5, D_sz=Dsz, npts=int(w.sum())))
print(f"crossmatched with usable window: {len(rows)} clusters ({skipped} skipped <3 pts)")
assert len(rows) >= 10, "overlap too small"

D = np.array([r["D"] for r in rows]); S = np.array([r["sig"] for r in rows]); M = np.array([r["M"] for r in rows])
Dsz = np.array([r["D_sz"] for r in rows if r["D_sz"] is not None])
print(f"\n--- the mixed-method deviation (caustic g_obs x X-ray g_bar) ---")
for lab, m in (("sigma > 1000", S > 1000), ("700-1000", (S > 700) & (S <= 1000)), ("<= 700", S <= 700)):
    if m.sum():
        print(f"  {lab:12s}: N={m.sum():3d}  D = {D[m].mean():+.3f} +- {D[m].std()/np.sqrt(m.sum()):.3f}")
print(f"  ALL         : N={len(D):3d}  D = {D.mean():+.3f} +- {D.std()/np.sqrt(len(D)):.3f}")
print(f"  forks: SZ-mass D = {Dsz.mean():+.3f} (N={len(Dsz)}); c=5 fork D = {np.mean([r['D_c5'] for r in rows]):+.3f}")
print(f"  frozen P2 threshold: normal-phase excess > +0.15;  T3 X-ray value: +0.44 (V06)")

# mass-at-fixed-sigma: partial Spearman of D vs log M controlling log sigma (and vice versa)
from scipy.stats import spearmanr
def partial(a, b, ctrl):
    ra = a - np.poly1d(np.polyfit(ctrl, a, 1))(ctrl)
    rb = b - np.poly1d(np.polyfit(ctrl, b, 1))(ctrl)
    return spearmanr(ra, rb)
pm = partial(D, np.log10(M), np.log10(S))
ps = partial(D, np.log10(S), np.log10(M))
print(f"\n--- the de-entangled mass-at-fixed-sigma test (deferred leg of P4) ---")
print(f"  partial rho(D, logM | logSigma) = {pm.statistic:+.2f} (p = {pm.pvalue:.3f})")
print(f"  partial rho(D, logSigma | logM) = {ps.statistic:+.2f} (p = {ps.pvalue:.3f})")
print(f"  phase-boundary reading expects: deviation organized by sigma (phase), flat in M at fixed sigma.")

verdict = dict(
    n=len(rows), D_all=float(D.mean()), D_err=float(D.std()/np.sqrt(len(D))),
    D_sz_fork=float(Dsz.mean()), D_c5_fork=float(np.mean([r['D_c5'] for r in rows])),
    partial_M_given_sigma=[float(pm.statistic), float(pm.pvalue)],
    partial_sigma_given_M=[float(ps.statistic), float(ps.pvalue)],
    p2_mixed="PASS" if D.mean() > 0.15 else "FAIL",
    note="g_obs caustic (galaxy dynamics), g_bar ACCEPT entropy+<kT> reconstruction + nominal BCG; "
         "isothermal-kT approximation and mass-per-electron 1.97e-24 g logged as conventions")
json.dump(dict(rows=rows, verdict=verdict), open(os.path.join(DERIVED, "t3_7_mixed_method.json"), "w"), indent=1)
print(f"\nP2 on the mixed-method sample: {verdict['p2_mixed']} (D = {D.mean():+.3f})")
print("wrote derived/t3_7_mixed_method.json")
