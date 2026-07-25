#!/usr/bin/env python
"""
T3.5 — systematics: refit sigma_b (and the m inversion) under the named forks.
 V0 baseline (reproduces T3.3); V1 maximal-baryon fork (hot-halo systems' g_bar ->
 0.156*g_obs — upper bound on feedback contamination); V2 cross-method offset
 (S09 D shifted by the measured overlap offset -0.06, worst-case -0.13);
 V3 stellar-budget forks (aperture x M/L: max-stars 1.0x1.2, min-stars 0.5x0.7);
 V4 quality cut (drop tier-2 and estimated-T groups).
Combined budget: m = fiducial-inversion(sigma_b) with stat CI + syst envelope + rho fork.
"""
import os, sys, re, json, io, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from t3_guard import enforce, ROOT, DERIVED

freezes, _ = enforce(staged=[os.path.join(ROOT, "data/T3/sun09/ms.tex"),
                             os.path.join(ROOT, "data/T3/vikhlinin06/mprof.tex")],
                     allow_dirty_paths=("t3_5_systematics.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
FB = 0.156
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19; MP = 1.6726e-27
rho_c0 = 3*(H0*1e3/(KPC*1e3))**2/(8*np.pi*G)
g_rar = lambda gb: gb/(1 - np.exp(-np.sqrt(gb/A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\\tablenotemark\{[a-z]\}|~", " ",
                         re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")

# ---------- S09: per-point (gobs, Mgas_frac, Mstar, M) so forks recompute ----------
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
    T500 = num1(cells[1]); est_T = "(" in cells[1]; tier2 = "*" in cells[1]
    r500 = num1(cells[3]); r2500 = num1(cells[4]); M500t = num1(cells[5])
    fg500 = num1(cells[6]); fg2500 = num1(cells[7])
    if np.isnan(T500) or np.isnan(r2500) or np.isnan(fg2500): continue
    rho_cz = rho_c0*(OM*(1+z)**3 + 1-OM)
    Ms = (10**LK)*MSUN if not np.isnan(LK) else 0.0
    apts = []
    r_m = r2500*KPC; M = (4*np.pi/3)*2500*rho_cz*r_m**3
    apts.append(dict(gobs=G*M/r_m**2, fgas=fg2500, Ms=Ms, M=M, r=r_m, ap=0.7))
    if not np.isnan(fg500) and not np.isnan(r500):
        r_m = r500*KPC
        M = M500t*1e13*MSUN if not np.isnan(M500t) else (4*np.pi/3)*500*rho_cz*r_m**3
        apts.append(dict(gobs=G*M/r_m**2, fgas=fg500, Ms=Ms, M=M, r=r_m, ap=1.0))
    S09.append(dict(name=name, sigma=400.4*np.sqrt(T500), est_T=est_T, tier2=tier2, apts=apts))

def s09_D(g, fork_fb=False, ap_f=1.0, ml=1.0, offset=0.0):
    pts = []
    for a in g["apts"]:
        gbar = a["gobs"]*FB if fork_fb else G*(a["fgas"]*a["M"] + ap_f*a["ap"]*ml*a["Ms"])/a["r"]**2
        if GLO <= gbar <= GHI:
            pts.append(np.log10(a["gobs"]) - np.log10(g_rar(gbar)))
    return (np.mean(pts) + offset) if pts else None

# ---------- V06: profiles (rg, gobs, gbar) for refork ----------
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
    dm = t_d[name][1:10]
    vals = t_m[name] + [np.nan]*9
    r500, Tspec, c500, M2500 = vals[0], vals[1], vals[3], vals[4]
    if any(np.isnan(x) for x in (r500, c500, M2500, Tspec)): continue
    rho_cz = rho_c0*(OM*(1+z)**3 + 1-OM)
    r2500 = (3*M2500*1e14*MSUN/(4*np.pi*2500*rho_cz))**(1/3)/KPC
    rs_n = r500/c500; mu = lambda x: np.log(1+x)-x/(1+x)
    rg = np.geomspace(max(rmin, 8), rdet, 200)
    Mn = M2500*1e14*MSUN*mu(rg/rs_n)/mu(r2500/rs_n)
    rho_g = 1.624*MP*np.sqrt(npne(rg, dm)*1e12)
    Mg = np.concatenate([[0], np.cumsum(0.5*(rho_g[1:]*(rg[1:]*KPC)**2 + rho_g[:-1]*(rg[:-1]*KPC)**2)*4*np.pi*np.diff(rg)*KPC)])
    gobs = G*Mn/(rg*KPC)**2; gbar = G*(Mg + 1e12*MSUN)/(rg*KPC)**2
    V06.append(dict(name=name, sigma=400.4*np.sqrt(Tspec), gobs=gobs, gbar=gbar))
def v06_D(c, fork_fb=False):
    gb = c["gobs"]*FB if fork_fb else c["gbar"]
    m = (gb >= GLO) & (gb <= GHI)
    return float(np.mean(np.log10(c["gobs"][m]) - np.log10(g_rar(gb[m])))) if m.sum() >= 3 else None

# ---------- SPARC + field ETGs (fixed classes) ----------
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
def run(tag, s09_kw=None, v06_fork=False, s09_filter=None):
    pts = list(base_pts)
    for g in S09:
        if s09_filter and not s09_filter(g): continue
        d = s09_D(g, **(s09_kw or {}))
        if d is not None: pts.append((g["sigma"], d, "Sun09"))
    for c in V06:
        d = v06_D(c, fork_fb=v06_fork)
        if d is not None: pts.append((c["sigma"], d, "V06"))
    sig = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
    w = np.array([CLS_SIG[p[2]] for p in pts])
    p, _ = curve_fit(model, np.log(sig), D, p0=[0., .45, np.log(400.), .3],
                     sigma=w, bounds=([-.2, 0, np.log(100), .02], [.2, 1, np.log(2000), 1.5]), maxfev=30000)
    sb = float(np.exp(p[2]))
    m = FZ["kernel"]["m_eV"]["fiducial"]*(sb/600.)**(-0.75)
    print(f"  {tag:44s} sigma_b = {sb:6.0f} km/s   m = {m:5.1f} eV   (step {p[1]:.2f}, width {p[3]:.2f})")
    return dict(tag=tag, sigma_b=sb, m_eV=float(m), step=float(p[1]), width=float(p[3]))

print("=== T3.5 SYSTEMATICS FOREST ===")
res = []
res.append(run("V0 baseline"))
res.append(run("V1 maximal-baryon fork (S09+V06 hot halos)", s09_kw=dict(fork_fb=True), v06_fork=True))
res.append(run("V2 method offset: S09 D - 0.06", s09_kw=dict(offset=-0.06)))
res.append(run("V2b method offset worst-case: S09 D - 0.13", s09_kw=dict(offset=-0.13)))
res.append(run("V3 max-stars (aperture 1.0, M/L 1.2)", s09_kw=dict(ap_f=1.0/0.7, ml=1.2)))
res.append(run("V3b min-stars (aperture 0.5, M/L 0.7)", s09_kw=dict(ap_f=0.5/0.7, ml=0.7)))
res.append(run("V4 quality: drop tier-2 & estimated-T", s09_filter=lambda g: not (g["tier2"] or g["est_T"])))

sbs = [r["sigma_b"] for r in res]; ms = [r["m_eV"] for r in res]
stat = json.load(open(os.path.join(DERIVED, "t3_3_straddle.json")))["sigma_b"]
rho_f = FZ["kernel"]["rho_convention_fork"]
m_stat = (FZ["kernel"]["m_eV"]["fiducial"]*(stat["hi84"]/600.)**(-0.75),
          FZ["kernel"]["m_eV"]["fiducial"]*(stat["lo16"]/600.)**(-0.75))
m_sys = (min(ms), max(ms))
m_out = (m_stat[0]*rho_f[0]**0.25*min(1, m_sys[0]/res[0]["m_eV"]),
         m_stat[1]*rho_f[1]**0.25*max(1, m_sys[1]/res[0]["m_eV"]))
print("\n=== COMBINED BUDGET ===")
print(f"sigma_b: central {res[0]['sigma_b']:.0f}, stat [{stat['lo16']:.0f}, {stat['hi84']:.0f}], syst envelope [{min(sbs):.0f}, {max(sbs):.0f}]")
print(f"m: central {res[0]['m_eV']:.1f} eV, stat [{m_stat[0]:.1f}, {m_stat[1]:.1f}], syst x[{m_sys[0]/res[0]['m_eV']:.2f}, {m_sys[1]/res[0]['m_eV']:.2f}], rho x[{rho_f[0]**0.25:.2f}, {rho_f[1]**0.25:.2f}]")
print(f"FINAL: m = {res[0]['m_eV']:.0f} eV, combined band [{m_out[0]:.0f}, {m_out[1]:.0f}] eV  vs frozen window [{FZ['kernel']['m_eV']['min']:.1f}, {FZ['kernel']['m_eV']['max']:.1f}]")
p2f = np.mean([v06_D(c, fork_fb=True) for c in V06 if v06_D(c, fork_fb=True) is not None])
print(f"P2 under maximal fork (V06 clusters): D' = {p2f:+.3f}  ({'still PASS' if p2f > 0.15 else 'FAILS'})")

json.dump(dict(frozen_against=FZ["provenance"], variants=res,
               stat_sigma_b=stat, m_final=dict(central=res[0]["m_eV"], band=[m_out[0], m_out[1]]),
               P2_under_fork=float(p2f)),
          open(os.path.join(DERIVED, "t3_5_systematics.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.2, 3.2))
for i, r_ in enumerate(res):
    y = len(res)-i
    ax.plot(r_["sigma_b"], y, "o", ms=6, color="crimson" if i == 0 else ".3")
    ax.text(120, y, r_["tag"], fontsize=7, va="center", ha="right")
    ax.text(r_["sigma_b"]*1.06, y, f"m={r_['m_eV']:.1f}", fontsize=6.5, va="center")
ax.axvspan(stat["lo16"], stat["hi84"], color="crimson", alpha=.10, lw=0)
ax.axvspan(476, 756, color="orange", alpha=.10, lw=0)
ax.axvline(600, color="orange", ls=":", lw=1)
ax.set(xscale="log", xlim=(100, 1000), ylim=(.3, len(res)+.7), yticks=[],
       xlabel=r"$\sigma_b$ [km s$^{-1}$]")
ax.set_title("T3.5: break location under the systematics forks", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_t3_forks.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/t3_5_systematics.json, figures/fig_t3_forks.png")