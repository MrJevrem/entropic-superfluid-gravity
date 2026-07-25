#!/usr/bin/env python
"""
T3.6 (slot reassigned from skipped GAMA item) — the f_b(sigma) partial-restoration model.
The principled middle between T3.5's V0 (no correction) and V1 (cosmic-f_b, over-generous:
even retentive clusters don't reach cosmic f_b at r2500 because gas is less concentrated
than DM). Baseline = the RETENTIVE value: mean f_b at each aperture over the hottest
systems (sigma > 900, both samples), which feedback cannot have depleted. Correction:
f_eff = max(f_meas, f_base) per aperture — lift only below-baseline systems, never
touch clusters. Refit sigma_b -> the pinned midpoint and the tightened m band.
Known one-sided bias flagged: field-ETG budgets (stars-only) omit hot-gas atmospheres
-> shelf D biased HIGH -> true sigma_b biased LOW by that channel (direction stated).
"""
import os, sys, re, json, io, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from t3_guard import enforce, ROOT, DERIVED

freezes, _ = enforce(staged=[os.path.join(ROOT, "data/T3/sun09/ms.tex"),
                             os.path.join(ROOT, "data/T3/vikhlinin06/mprof.tex")],
                     allow_dirty_paths=("t3_6_fb_baseline.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19; MP = 1.6726e-27
rho_c0 = 3*(H0*1e3/(KPC*1e3))**2/(8*np.pi*G)
g_rar = lambda gb: gb/(1 - np.exp(-np.sqrt(gb/A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\\tablenotemark\{[a-z]\}|~", " ",
                         re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")

# ---- S09 (as t3_5: per-aperture gobs, f_meas) ----
tex = open(os.path.join(ROOT, "data/T3/sun09/ms.tex")).read()
zLK = {}
for line in tex.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) >= 8 and re.search(r"\d", cells[1]):
        name = norm(cells[0]); z, LK = num1(cells[1]), num1(cells[7])
        if name and not np.isnan(z): zLK[name] = (z, LK)
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
    Ms = (10**LK)*MSUN if not np.isnan(LK) else 0.0
    apts = []
    r_m = r2500*KPC; M = (4*np.pi/3)*2500*rho_cz*r_m**3
    apts.append(dict(gobs=G*M/r_m**2, f=fg2500 + 0.7*Ms/M, ap="2500"))
    if not np.isnan(fg500) and not np.isnan(r500):
        r_m = r500*KPC
        M5 = M500t*1e13*MSUN if not np.isnan(M500t) else (4*np.pi/3)*500*rho_cz*r_m**3
        apts.append(dict(gobs=G*M5/r_m**2, f=fg500 + Ms/M5, ap="500"))
    S09.append(dict(name=name, sigma=400.4*np.sqrt(T500), apts=apts))

# ---- V06 (fb at 2500 from their table; profile D scaled by baseline lift) ----
vtex = open(os.path.join(ROOT, "data/T3/vikhlinin06/mprof.tex")).read().replace(r"\pz", " ").replace("~", " ")
blocks = re.findall(r"\\startdata(.*?)\\enddata", vtex, re.S)
def vrows(block):
    out = []
    for line in re.split(r"\\\\", block):
        line = re.sub(r"\\mcc\{[^}]*\}", "nan", line.replace("\\dotfill", "").replace("$", "")).replace(r"\nodata", "nan")
        cells = [c.strip() for c in re.sub(r"\\,", "", line).split("&")]
        if len(cells) < 3: continue
        name = re.sub(r"[\s{}\\]", "", cells[0]); vals = []
        for c in cells[1:]:
            mm = re.match(r"^\s*([-+]?\d*\.?\d+)", c)
            vals.append(float(mm.group(1)) if mm else np.nan)
        out.append((name, vals))
    return out
t_s = {n: v for n, v in vrows(blocks[0])}; t_d = {n: v for n, v in vrows(blocks[1])}; t_m = {n: v for n, v in vrows(blocks[2])}
def npne(r, p):
    n0, rc, rs, al, be, ep, n02, rc2, b2 = p
    t1 = (n0*1e-3)**2*(r/rc)**(-al)/(1+r**2/rc**2)**(3*be-al/2)/(1+r**3/rs**3)**(ep/3)
    t2 = 0.0 if np.isnan(n02) else (n02*1e-1)**2/(1+r**2/rc2**2)**(3*b2)
    return t1 + t2
V06 = []
for name, sv in t_s.items():
    if name not in t_d or name not in t_m: continue
    z, rmin, rdet = sv[0], sv[1], sv[2]
    vals = t_m[name] + [np.nan]*9
    r500, Tspec, c500, M2500, fg2500 = vals[0], vals[1], vals[3], vals[4], vals[6]
    if any(np.isnan(x) for x in (r500, c500, M2500, Tspec)): continue
    rho_cz = rho_c0*(OM*(1+z)**3 + 1-OM)
    r2500 = (3*M2500*1e14*MSUN/(4*np.pi*2500*rho_cz))**(1/3)/KPC
    rs_n = r500/c500; mu = lambda x: np.log(1+x)-x/(1+x)
    rg = np.geomspace(max(rmin, 8), rdet, 200)
    Mn = M2500*1e14*MSUN*mu(rg/rs_n)/mu(r2500/rs_n)
    rho_g = 1.624*MP*np.sqrt(npne(rg, t_d[name][1:10])*1e12)
    Mg = np.concatenate([[0], np.cumsum(0.5*(rho_g[1:]*(rg[1:]*KPC)**2 + rho_g[:-1]*(rg[:-1]*KPC)**2)*4*np.pi*np.diff(rg)*KPC)])
    gobs = G*Mn/(rg*KPC)**2; gbar = G*(Mg + 1e12*MSUN)/(rg*KPC)**2
    fb2500 = (fg2500 + 0.7*1e12*MSUN/(M2500*1e14*MSUN)) if not np.isnan(fg2500) else np.nan
    V06.append(dict(name=name, sigma=400.4*np.sqrt(Tspec), gobs=gobs, gbar=gbar, fb2500=fb2500))

# ---- the baseline from the hot set ----
hot_f25 = [a["f"] for g in S09 if g["sigma"] > 900 for a in g["apts"] if a["ap"] == "2500"] + \
          [c["fb2500"] for c in V06 if c["sigma"] > 900 and not np.isnan(c["fb2500"])]
hot_f50 = [a["f"] for g in S09 if g["sigma"] > 900 for a in g["apts"] if a["ap"] == "500"]
FB25 = float(np.mean(hot_f25)); FB25_lo = FB25 - float(np.std(hot_f25))
FB50 = float(np.mean(hot_f50)) if hot_f50 else 0.115
print(f"retentive baselines (sigma>900): f_b(2500) = {FB25:.3f} +/- {np.std(hot_f25):.3f} (N={len(hot_f25)}); "
      f"f_b(500) = {FB50:.3f} (N={len(hot_f50)})")
print(f"vs cosmic 0.156: even retentive clusters reach only {FB25/0.156:.0%} at r2500 — V1's overcorrection quantified")

def D_s09(g, fb25, fb50):
    pts = []
    for a in g["apts"]:
        base = fb25 if a["ap"] == "2500" else fb50
        gbar = a["gobs"]*max(a["f"], base)
        if GLO <= gbar <= GHI: pts.append(np.log10(a["gobs"]) - np.log10(g_rar(gbar)))
    return np.mean(pts) if pts else None
def D_v06(c, fb25):
    lift = max(1.0, fb25/c["fb2500"]) if not np.isnan(c["fb2500"]) else 1.0
    gb = c["gbar"]*lift
    m = (gb >= GLO) & (gb <= GHI)
    return float(np.mean(np.log10(c["gobs"][m]) - np.log10(g_rar(gb[m])))) if m.sum() >= 3 else None

# ---- assemble + fit (as t3_5) ----
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
    m = (gbar > 0) & (gbar >= GLO) & (gbar <= GHI)
    if m.sum() >= 3:
        base_pts.append((float(np.mean(vobs[-3:]))/np.sqrt(2),
                         float(np.mean(np.log10(gobs[m]) - np.log10(g_rar(gbar[m])))), "SPARC"))
for r_ in json.load(open(os.path.join(DERIVED, "t3_2_shelf.json")))["rows"]:
    if r_["rho_env"] < 0.5: base_pts.append((r_["sigma"], r_["D"], "ETG_field"))
CLS_SIG = {"SPARC": .20, "ETG_field": .22, "Sun09": .12, "V06": .08}
def model(lnS, a, b, lnSb, wid): return a + b/(1+np.exp(-(lnS-lnSb)/wid))
def fit_once(fb25, fb50, boot_rng=None):
    pts = list(base_pts)
    for g in S09:
        d = D_s09(g, fb25, fb50)
        if d is not None: pts.append((g["sigma"], d, "Sun09"))
    for c in V06:
        d = D_v06(c, fb25)
        if d is not None: pts.append((c["sigma"], d, "V06"))
    sig = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
    cls = np.array([p[2] for p in pts]); w = np.array([CLS_SIG[c] for c in cls])
    idx = np.arange(len(sig))
    if boot_rng is not None:
        idx = np.concatenate([boot_rng.choice(np.where(cls == c)[0], (cls == c).sum(), replace=True)
                              for c in CLS_SIG])
    p, _ = curve_fit(model, np.log(sig[idx]), D[idx], p0=[0., .4, np.log(400.), .3],
                     sigma=w[idx], bounds=([-.2, 0, np.log(100), np.log(2000)*0+.02], [.2, 1, np.log(2000), 1.5]),
                     maxfev=30000)
    return p, pts

M_FID = FZ["kernel"]["m_eV"]["fiducial"]
m_of = lambda s: M_FID*(s/600.0)**(-0.75)
print("\n=== T3.6: the pinned midpoint ===")
res = {}
for tag, f25, f50 in (("V0.5 baseline-mean", FB25, FB50), ("V0.5lo baseline-1sd", FB25_lo, FB50-0.02)):
    p, pts = fit_once(f25, f50)
    sb = float(np.exp(p[2]))
    res[tag] = dict(f25=f25, f50=f50, sigma_b=sb, m=float(m_of(sb)), step=float(p[1]), width=float(p[3]))
    print(f"  {tag:22s} f25={f25:.3f}: sigma_b = {sb:6.0f} km/s  m = {m_of(sb):5.1f} eV  (step {p[1]:.2f}, width {p[3]:.2f})")
rng = np.random.default_rng(20260726)
boots = []
for _ in range(200):
    try: boots.append(np.exp(fit_once(FB25, FB50, rng)[0][2]))
    except Exception: pass
sb_lo, sb_med, sb_hi = np.percentile(boots, [16, 50, 84])
print(f"  bootstrap (V0.5): sigma_b = {sb_med:.0f} [{sb_lo:.0f}, {sb_hi:.0f}] -> m = {m_of(sb_med):.1f} [{m_of(sb_hi):.1f}, {m_of(sb_lo):.1f}] eV")
t35 = json.load(open(os.path.join(DERIVED, "t3_5_systematics.json")))
print(f"\ncontext: V0 254 (m 16.1) | V0.5 {res['V0.5 baseline-mean']['sigma_b']:.0f} "
      f"(m {res['V0.5 baseline-mean']['m']:.1f}) | V1 829 (m 6.6, retained as over-generous bound)")
rho_f = FZ["kernel"]["rho_convention_fork"]
mlo = m_of(sb_hi)*rho_f[0]**0.25; mhi = m_of(sb_lo)*rho_f[1]**0.25
sys_lo = min(res["V0.5lo baseline-1sd"]["m"], m_of(sb_hi)); sys_hi = max(16.1, m_of(sb_lo))  # V0 kept as no-correction upper edge
print(f"proposed tightened band: m = {m_of(sb_med):.0f} eV, stat x rho [{mlo:.0f}, {mhi:.0f}], "
      f"principled-systematics [{sys_lo:.0f}, {sys_hi:.0f}]  (was [5, 25])")
print("ETG-shelf caveat: stars-only budgets omit hot-gas atmospheres -> shelf D biased high ->")
print("sigma_b from this fit biased LOW by that channel (direction stated; gas data not in hand)")

json.dump(dict(frozen_against=FZ["provenance"], baselines=dict(fb2500=FB25, fb2500_sd=float(np.std(hot_f25)),
               fb500=FB50, cosmic=0.156, retention_vs_cosmic=FB25/0.156),
               fits=res, bootstrap=dict(median=float(sb_med), lo=float(sb_lo), hi=float(sb_hi)),
               m_band_proposal=dict(central=float(m_of(sb_med)), stat_rho=[float(mlo), float(mhi)],
                                     syst=[float(sys_lo), float(sys_hi)], previous=[5, 25]),
               etg_caveat="stars-only shelf budgets omit hot gas: sigma_b biased low (direction stated)"),
          open(os.path.join(DERIVED, "t3_6_fb_baseline.json"), "w"), indent=1)
print("\nwrote derived/t3_6_fb_baseline.json")