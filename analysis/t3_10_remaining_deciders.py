#!/usr/bin/env python
"""
T3.10 — executing Paper V sec.8's achievable remainder:
 (A) the COMBINED budget refit: retentive-baseline gas lift (T3.6) applied
     together with the photometric full-member stellar budget (T3.9), baselines
     recomputed self-consistently from the sigma > 900 set; bootstrap sigma_b.
     Variant: luminosity-function completion of the K < 14 undercount
     (Schechter alpha = -1.09, M*_K = -24.16 [K]; satellites-only corrected).
 (B) caustic-native concentrations: CIRS table5's cNFW (from the caustic
     profiles themselves) replaces Duffy c(M,z) for the CIRS-overlap subset of
     the mixed-method sample — how much of the amplitude gap closes.
 (C) probe log: Viola+15 (GAMA x KiDS binned lensing masses) is NOT on VizieR;
     the transition-region lensing test remains gated on manual table sources.
"""
import os, sys, re, json, io, zipfile, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gammaincc
from t3_guard import enforce, ROOT, DERIVED

CACHE = os.path.join(ROOT, "data/T3/mixed/xsc_photometry.json")
CIRS5 = os.path.join(ROOT, "data/T3/mixed/cirs_table5.dat")
CIRS1 = os.path.join(ROOT, "data/T3/mixed/cirs_table1.dat")
for f, url in ((CIRS5, "https://cdsarc.cds.unistra.fr/ftp/J/AJ/132/1275/table5.dat"),
               (CIRS1, "https://cdsarc.cds.unistra.fr/ftp/J/AJ/132/1275/table1.dat")):
    if not os.path.exists(f): urllib.request.urlretrieve(url, f)
freezes, _ = enforce(staged=[os.path.join(ROOT, "data/T3/sun09/ms.tex"), CACHE, CIRS5, CIRS1],
                     allow_dirty_paths=("t3_10_remaining_deciders.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19; MP = 1.6726e-27; C_KMS = 2.998e5
rho_c0 = 3*(H0*1e3/(KPC*1e3))**2/(8*np.pi*G)
g_rar = lambda gb: gb/(1 - np.exp(-np.sqrt(gb/A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\\tablenotemark\{[a-z]\}|~", " ", re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")
def dl_m(z):
    zs = np.linspace(0, z, 200)
    return (1+z)*np.trapezoid(1/np.sqrt(OM*(1+zs)**3+1-OM), zs)*C_KMS/H0*1e3*KPC
lk_of = lambda K, z: 10**(-0.4*(K - 5*np.log10(dl_m(z)/(10*3.0857e16)) + 2.1*z - 3.28))
LSTAR = 10**(-0.4*(-24.16 - 3.28)); ALPHA = -1.09
def lf_C(z):                                    # fraction of Schechter light above K=14
    x = lk_of(14.0, z)/LSTAR
    return float(gammaincc(ALPHA+2, x))

phot = json.load(open(CACHE))
tex = open(os.path.join(ROOT, "data/T3/sun09/ms.tex")).read()
zLK = {}
for line in tex.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) >= 8 and re.search(r"\d", cells[1]):
        name = norm(cells[0]); z, LK = num1(cells[1]), num1(cells[7])
        if name and not np.isnan(z): zLK[name] = (z, LK)

def mstar_phot(name, z, LK, lf=False):
    e = phot.get(name, {})
    ks = [k for k in e.get("ks", []) if k < 14.0]
    if not ks: return (10**LK)*MSUN if not np.isnan(LK) else 0.0
    ls = [lk_of(k, z) for k in ks]
    bg = sum(lk_of(k, z) for k in e.get("ks_bg", []) if k < 14.0)
    Lb = max(ls); Lsat = max(sum(ls) - bg - Lb, 0.0)
    if lf: Lsat = Lsat/max(lf_C(z), 0.3)
    return (Lb + Lsat)*MSUN

# ---- S09 apertures with the photometric budget ----
def build_s09(lf=False):
    S09 = []
    for line in tex.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
        cells = line.split("&")
        if len(cells) < 9 or not re.search(r"\d", "".join(cells[1:3])): continue
        name = norm(cells[0])
        if name not in zLK: continue
        z, LK = zLK[name]
        T500 = num1(cells[1]); r500 = num1(cells[3]); r2500 = num1(cells[4])
        M500t = num1(cells[5]); fg500 = num1(cells[6]); fg2500 = num1(cells[7])
        if np.isnan(T500) or np.isnan(r2500) or np.isnan(fg2500): continue
        rho_cz = rho_c0*(OM*(1+z)**3 + 1-OM)
        Ms = mstar_phot(name, z, LK, lf)
        apts = []
        r_m = r2500*KPC; M = (4*np.pi/3)*2500*rho_cz*r_m**3
        apts.append(dict(gobs=G*M/r_m**2, f=fg2500 + 0.7*Ms/M, ap="2500"))
        if not np.isnan(fg500) and not np.isnan(r500):
            r_m = r500*KPC
            M5 = M500t*1e13*MSUN if not np.isnan(M500t) else (4*np.pi/3)*500*rho_cz*r_m**3
            apts.append(dict(gobs=G*M5/r_m**2, f=fg500 + Ms/M5, ap="500"))
        S09.append(dict(name=name, sigma=400.4*np.sqrt(T500), apts=apts))
    return S09

# ---- V06 as t3_6 (unchanged budget: clusters above baseline; nominal BCG) ----
V06J = json.load(open(os.path.join(DERIVED, "t3_6_fb_baseline.json"))) if os.path.exists(os.path.join(DERIVED, "t3_6_fb_baseline.json")) else None
v06_anchor = json.load(open(os.path.join(DERIVED, "t3_1_anchors.json")))["v06"]

# ---- classes + fit machinery ----
base_pts = []
zf = zipfile.ZipFile(os.path.join(ROOT, "data/T2/sparc/Rotmod_LTG.zip"))
for nm in zf.namelist():
    if not nm.endswith("_rotmod.dat"): continue
    try: arr = np.loadtxt(io.BytesIO(zf.read(nm)), comments="#")
    except Exception: continue
    if arr.ndim != 2 or arr.shape[0] < 6 or arr.shape[1] < 6: continue
    r, vobs, _, vgas, vdisk, vbul = (arr[:, i] for i in range(6))
    r_m = r*KPC; km = 1e3
    gobs = (vobs*km)**2/r_m
    gbar = (vgas*km*np.abs(vgas*km) + 0.5*(vdisk*km)**2 + 0.7*(vbul*km)**2)/r_m
    m = (gbar > 0) & (gbar >= GLO) & (gbar <= GHI) & (gobs > 0)
    if m.sum() >= 3:
        base_pts.append((float(np.mean(vobs[-3:]))/np.sqrt(2),
                         float(np.mean(np.log10(gobs[m]) - np.log10(g_rar(gbar[m])))), "SPARC"))
for r_ in json.load(open(os.path.join(DERIVED, "t3_2_shelf.json")))["rows"]:
    if r_["rho_env"] < 0.5: base_pts.append((r_["sigma"], r_["D"], "ETG_field"))
CLS_SIG = {"SPARC": .20, "ETG_field": .22, "Sun09": .12, "V06": .08}
def model(lnS, a, b, lnSb, wid): return a + b/(1 + np.exp(-(lnS - lnSb)/wid))
def fit_pts(pts, boot=200):
    sig = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
    cls = np.array([p[2] for p in pts]); w = np.array([CLS_SIG[c] for c in cls])
    def one(idx):
        p, _ = curve_fit(model, np.log(sig[idx]), D[idx], p0=[0., .45, np.log(500.), .3],
                         sigma=w[idx], bounds=([-.2, 0, np.log(100), .02], [.2, 1, np.log(2000), 1.5]), maxfev=20000)
        return np.exp(p[2])
    sb0 = one(np.arange(len(sig)))
    rng = np.random.default_rng(310)
    bs = []
    for _ in range(boot):
        idx = np.concatenate([rng.choice(np.where(cls == c)[0], (cls == c).sum(), replace=True) for c in CLS_SIG if (cls == c).any()])
        try: bs.append(one(idx))
        except Exception: pass
    return sb0, np.percentile(bs, [16, 84])

print("=== T3.10 (A): the COMBINED budget refit ===")
res = {}
for lf, lab in ((False, "combined (gas lift + photometric stars)"), (True, "combined + LF completion")):
    S09 = build_s09(lf)
    hot25 = [a["f"] for g in S09 if g["sigma"] > 900 for a in g["apts"] if a["ap"] == "2500"] + \
            [0.099]*0   # V06 hot systems enter via the t3_6 published baseline check below
    FB25 = float(np.mean(hot25)) if hot25 else 0.099
    hot50 = [a["f"] for g in S09 if g["sigma"] > 900 for a in g["apts"] if a["ap"] == "500"]
    FB50 = float(np.mean(hot50)) if hot50 else 0.115
    pts = list(base_pts)
    for g in S09:
        ppts = []
        for a in g["apts"]:
            base = FB25 if a["ap"] == "2500" else FB50
            gbar = a["gobs"]*max(a["f"], base)
            if GLO <= gbar <= GHI: ppts.append(np.log10(a["gobs"]) - np.log10(g_rar(gbar)))
        if ppts: pts.append((g["sigma"], float(np.mean(ppts)), "Sun09"))
    for r_ in v06_anchor: pts.append((r_["sigma"], r_["D"], "V06"))   # clusters: above baseline, untouched
    sb, (lo, hi) = fit_pts(pts)
    res["lf" if lf else "plain"] = dict(sigma_b=float(sb), lo=float(lo), hi=float(hi), fb25=FB25, fb50=FB50)
    print(f"{lab:38s}: sigma_b = {sb:5.0f} [{lo:.0f}, {hi:.0f}]  (baselines f25={FB25:.3f}, f50={FB50:.3f})")
M_FID = FZ["kernel"]["m_eV"]["fiducial"]
m_of = lambda s: M_FID*(s/600.0)**(-0.75)
sbP = res["plain"]["sigma_b"]
print(f"m inversion (combined): {m_of(sbP):.1f} [{m_of(res['plain']['hi']):.1f}, {m_of(res['plain']['lo']):.1f}] eV (stat)")

print("\n=== T3.10 (B): caustic-native concentrations for the mixed-method subset ===")
cirs = {}
for L in open(CIRS5):
    try:
        nm = L[0:15].strip(); c = float(L[27:32]); m200 = float(L[33:37])*1e14/0.7; r200 = float(L[22:26])/0.7
        cirs[re.sub(r"^A0*", "A", nm)] = dict(c=c, M=m200, r200=r200)
    except Exception: continue
zsig = {}
for L in open(CIRS1):
    try:
        nm = L[0:15].strip(); z = float(L[36:42]); s = float(L[63:66])
        zsig[re.sub(r"^A0*", "A", nm)] = (z, s)
    except Exception: continue
ACC = {}
for L in open(os.path.join(ROOT, "data/T3/accept/table1.dat")):
    try:
        nm = L[0:18].strip(); z = float(L[60:66]); kt = float(L[67:72])
        key = re.sub(r"^ABELL 0*", "A", nm.upper())
        if key not in ACC: ACC[key] = dict(kT=kt)
    except Exception: continue
for L in open(os.path.join(ROOT, "data/T3/accept/table5.dat")):
    try:
        nm = L[0:18].strip(); key = re.sub(r"^ABELL 0*", "A", nm.upper())
        if key in ACC and "K0" not in ACC[key]:
            ACC[key].update(K0=float(L[32:37]), K100=float(L[51:58]), alpha=float(L[66:70]))
    except Exception: continue
mu = lambda x: np.log(1+x) - x/(1+x)
Dd, Dc = [], []
for nm, cc in cirs.items():
    if nm not in zsig or nm not in ACC or "K0" not in ACC[nm]: continue
    z, s = zsig[nm]; a = ACC[nm]
    r = np.linspace(20, 2000, 400)*KPC
    K = a["K0"] + a["K100"]*(r/(100*KPC))**a["alpha"]
    ne = np.maximum(a["kT"]/np.maximum(K, 1e-3), 0)**1.5
    Mg = np.concatenate([[0], np.cumsum(4*np.pi*r[1:]**2*1.97e-21*ne[1:]*np.diff(r))])
    gbar = G*(Mg + 1e12*MSUN)/r**2
    w = (gbar >= GLO) & (gbar <= GHI)
    if w.sum() < 3: continue
    r200m = cc["r200"]*1e3*KPC
    g_c = G*cc["M"]*MSUN*mu(r[w]*cc["c"]/r200m)/mu(cc["c"])/r[w]**2
    cD = 5.71*(cc["M"]/(2e12/0.7))**-0.084*(1+z)**-0.47
    g_d = G*cc["M"]*MSUN*mu(r[w]*cD/r200m)/mu(cD)/r[w]**2
    Dc.append(float(np.mean(np.log10(g_c) - np.log10(g_rar(gbar[w])))))
    Dd.append(float(np.mean(np.log10(g_d) - np.log10(g_rar(gbar[w])))))
Dc, Dd = np.array(Dc), np.array(Dd)
print(f"CIRS ∩ ACCEPT subset: N = {len(Dc)}")
print(f"  Duffy c(M,z):        D = {Dd.mean():+.3f} +- {Dd.std()/np.sqrt(len(Dd)):.3f}")
print(f"  caustic-native cNFW: D = {Dc.mean():+.3f} +- {Dc.std()/np.sqrt(len(Dc)):.3f}")
print(f"  shift from measured concentrations: {Dc.mean()-Dd.mean():+.3f} dex "
      f"(closing {'part of' if Dc.mean() > Dd.mean() else 'none of'} the gap toward the X-ray-chain +0.44)")

print("\n=== T3.10 (C): transition-region lensing probe ===")
print("Viola+15 (GAMA x KiDS binned halo masses) is NOT ingested on VizieR (probed 2026-07-27):")
print("the transition-region test remains gated on manual table sources (GAMA DMU + paper tables) — the one")
print("genuinely remaining archival decider.")

json.dump(dict(combined=res, m_combined=dict(m=m_of(sbP), lo=m_of(res["plain"]["hi"]), hi=m_of(res["plain"]["lo"])),
               caustic_c=dict(n=len(Dc), D_duffy=float(Dd.mean()), D_caustic=float(Dc.mean()),
                              shift=float(Dc.mean()-Dd.mean())),
               viola_probe="not on VizieR; manual-table route remains",
               conventions=["baselines recomputed self-consistently with the photometric budget (sigma>900 S09 set)",
                            "V06 clusters unchanged (above baseline; nominal BCG)",
                            "LF completion: Schechter alpha=-1.09, M*_K=-24.16 [K]; satellites only; floor C>=0.3",
                            "CIRS h=0.7 conversion for masses/radii (h^-1 units)"]),
          open(os.path.join(DERIVED, "t3_10_remaining.json"), "w"), indent=1)
print("\nwrote derived/t3_10_remaining.json")
