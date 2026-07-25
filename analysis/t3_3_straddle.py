#!/usr/bin/env python
"""
T3.3 — the straddle scan and break fit.
New data: Sun+09 (43 Chandra groups; staged e-print, parsed at runtime).
Per group: g_obs(Delta) = (4pi/3) Delta rho_c(z) G r_Delta  (pure geometry+z),
g_bar from their measured f_gas + their measured member K-band light (M/L_K=1,
matching the SLUGGS convention; aperture factors 0.7 (r2500) / 1.0 (r500), logged).
sigma = 400.4 sqrt(T500) (frozen host-system proxy).
Break fit: logistic in ln sigma over {SPARC, field ETGs (rho_env<0.5 mechanical rule),
Sun09, V06}; bootstrap -> sigma_b posterior -> m posterior via the frozen inversion.
"""
import os, sys, re, json, io, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from t3_guard import enforce, ROOT, DERIVED

TEX = os.path.join(ROOT, "data/T3/sun09/ms.tex")
freezes, _ = enforce(staged=[TEX], allow_dirty_paths=("t3_3_straddle.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
rho_c0 = 3 * (H0*1e3/(KPC*1e3))**2 / (8*np.pi*G)
g_rar = lambda gb: gb / (1 - np.exp(-np.sqrt(gb / A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\^\{[^}]*\}|_\{[^}]*\}|\\tablenotemark\{[a-z]\}|~", " ",
                         re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")

tex = open(TEX).read()
# --- Table 1: name -> z, L_Ks
zLK = {}
for line in tex.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) >= 8 and re.search(r"\d", cells[1]):
        name = norm(cells[0])
        z, LK = num1(cells[1]), num1(cells[7])
        if name and not np.isnan(z): zLK[name] = (z, LK)
# --- Table: derived properties I
rows09 = []
for line in tex.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) < 9 or not re.search(r"\d", "".join(cells[1:3])): continue
    name = norm(cells[0])
    if name not in zLK: continue
    z, LK = zLK[name]
    T500 = num1(cells[1]); est_T = "(" in cells[1]; tier2 = "*" in cells[1]
    r500 = num1(cells[3]); r2500 = num1(cells[4])
    M500tab = num1(cells[5]); fg500 = num1(cells[6]); fg2500 = num1(cells[7])
    if np.isnan(T500) or np.isnan(r2500) or np.isnan(fg2500): continue
    rho_cz = rho_c0 * (OM*(1+z)**3 + 1-OM)
    Mstar = (10**LK) * MSUN if not np.isnan(LK) else 0.0
    pts = []
    # Delta = 2500 (primary, measured)
    r_m = r2500 * KPC
    M2500 = (4*np.pi/3) * 2500 * rho_cz * r_m**3
    gobs = G * M2500 / r_m**2
    gbar = G * (fg2500 * M2500 + 0.7 * Mstar) / r_m**2
    if GLO <= gbar <= GHI: pts.append(np.log10(gobs) - np.log10(g_rar(gbar)))
    # Delta = 500 (where fgas500 measured)
    if not np.isnan(fg500) and not np.isnan(r500):
        r_m = r500 * KPC
        M500 = M500tab*1e13*MSUN if not np.isnan(M500tab) else (4*np.pi/3)*500*rho_cz*r_m**3
        gobs5 = G * M500 / r_m**2
        gbar5 = G * (fg500 * M500 + Mstar) / r_m**2
        if GLO <= gbar5 <= GHI: pts.append(np.log10(gobs5) - np.log10(g_rar(gbar5)))
    if pts:
        rows09.append(dict(name=name, z=z, sigma=round(400.4*np.sqrt(T500),0), T500=T500,
                           D=round(float(np.mean(pts)),3), n_ap=len(pts), est_T=est_T, tier2=tier2))
print(f"Sun09: {len(rows09)} groups with window points (of {len(zLK)} sample)")

# --- assemble classes ---
pts = []   # (sigma, D, class)
# SPARC (recompute, as T3.1)
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
    ok = gbar > 0
    m = ok & (gbar >= GLO) & (gbar <= GHI) & (gobs > 0)
    if m.sum() >= 3:
        D = float(np.mean(np.log10(gobs[m]) - np.log10(g_rar(gbar[m]))))
        pts.append((float(np.mean(vobs[-3:]))/np.sqrt(2), D, "SPARC"))
# field ETGs (mechanical rule rho_env < 0.5)
for r_ in json.load(open(os.path.join(DERIVED, "t3_2_shelf.json")))["rows"]:
    if r_["rho_env"] < 0.5:
        pts.append((r_["sigma"], r_["D"], "ETG_field"))
for r_ in rows09: pts.append((r_["sigma"], r_["D"], "Sun09"))
for r_ in json.load(open(os.path.join(DERIVED, "t3_1_anchors.json")))["v06"]:
    pts.append((r_["sigma"], r_["D"], "V06"))

CLS_SIG = {"SPARC": .20, "ETG_field": .22, "Sun09": .12, "V06": .08}
sig = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
cls = np.array([p[2] for p in pts]); w = np.array([CLS_SIG[c] for c in cls])

def model(lnS, a, b, lnSb, wid): return a + b/(1 + np.exp(-(lnS - lnSb)/wid))
def fit(idx):
    p, _ = curve_fit(model, np.log(sig[idx]), D[idx], p0=[0., .45, np.log(500.), .3],
                     sigma=w[idx], bounds=([-.2, 0, np.log(100), .02], [.2, 1, np.log(2000), 1.5]),
                     maxfev=20000)
    return p
p0 = fit(np.arange(len(sig)))
rng = np.random.default_rng(20260725)
boots = []
for _ in range(400):
    idx = np.concatenate([rng.choice(np.where(cls == c)[0], (cls == c).sum(), replace=True)
                          for c in CLS_SIG])
    try: boots.append(fit(idx))
    except Exception: pass
boots = np.array(boots)
sb = np.exp(boots[:, 2]); sb_med, sb_lo, sb_hi = np.median(sb), *np.percentile(sb, [16, 84])
M_FID = FZ["kernel"]["m_eV"]["fiducial"]
m_of = lambda s: M_FID * (s/600.0)**(-0.75)
m_med, m_lo, m_hi = m_of(sb_med), m_of(sb_hi), m_of(sb_lo)
rho_f = FZ["kernel"]["rho_convention_fork"]
m_outer = (m_lo*rho_f[0]**0.25, m_hi*rho_f[1]**0.25)
env = FZ["kernel"]["sigma_crit_km_s"]["outer_envelope"]

print("\n=== T3.3 BREAK FIT ===")
print(f"points: {len(sig)} ({', '.join(f'{c}:{(cls==c).sum()}' for c in CLS_SIG)})")
print(f"full fit: floor a = {p0[0]:+.3f}, step b = {p0[1]:.3f}, width = {p0[3]:.2f} (ln sigma)")
print(f"sigma_b = {sb_med:.0f} [{sb_lo:.0f}, {sb_hi:.0f}] km/s  (frozen fiducial 600, rho band 476-756)")
print(f"P3: inside outer envelope {env}: {'PASS' if env[0] < sb_med < env[1] else 'FAIL'}")
print(f"m posterior: {m_med:.1f} [{m_lo:.1f}, {m_hi:.1f}] eV (stat) x rho fork -> [{m_outer[0]:.1f}, {m_outer[1]:.1f}] eV")
print(f"(frozen window: [{FZ['kernel']['m_eV']['min']:.1f}, {FZ['kernel']['m_eV']['max']:.1f}] eV)")
t1 = np.array([not r_["tier2"] for r_ in rows09])
s09sig = np.array([r_["sigma"] for r_ in rows09]); s09D = np.array([r_["D"] for r_ in rows09])
print("\nSun09 D by sigma bin:")
for lo, hi in ((300, 400), (400, 500), (500, 600), (600, 700)):
    m = (s09sig >= lo) & (s09sig < hi)
    if m.any(): print(f"  {lo}-{hi}: N={m.sum():2d}  mean D = {s09D[m].mean():+.3f}  (scatter {s09D[m].std():.3f})")
print(f"robustness: excluding tier-2 S09 groups (N={int((~t1).sum())}) -> refit sigma_b = ", end="")
try:
    s09names = [r_["name"] for r_ in rows09]
    mask = np.ones(len(sig), bool)
    j = 0
    for i, c in enumerate(cls):
        if c == "Sun09":
            if rows09[j]["tier2"]: mask[i] = False
            j += 1
    p2 = fit(np.where(mask)[0]); print(f"{np.exp(p2[2]):.0f} km/s")
except Exception as e: print(f"(failed: {e})")

ov = {"A262": None, "MKW4": None, "A1991": None, "RXJ1159": None}
for r_ in rows09:
    if r_["name"] in ov: ov[r_["name"]] = r_["D"]
v06D = {x["name"].replace("+5531",""): x["D"] for x in json.load(open(os.path.join(DERIVED,"t3_1_anchors.json")))["v06"]}
print("\ncross-method overlaps (Sun09 fgas-aperture vs V06 profile):")
for k, v in ov.items():
    if v is not None and k in v06D: print(f"  {k:8s}: S09 D={v:+.3f}  V06 D={v06D[k]:+.3f}  delta={v-v06D[k]:+.3f}")

json.dump(dict(frozen_against=FZ["provenance"], n_points=len(sig),
               fit=dict(a=float(p0[0]), b=float(p0[1]), sigma_b=float(np.exp(p0[2])), width=float(p0[3])),
               sigma_b=dict(median=float(sb_med), lo16=float(sb_lo), hi84=float(sb_hi)),
               m_eV=dict(median=float(m_med), lo=float(m_lo), hi=float(m_hi),
                         rho_fork_outer=[float(m_outer[0]), float(m_outer[1])]),
               sun09=rows09,
               conventions=["M* = 10^LKs (M/L_K = 1, SLUGGS-matched); aperture 0.7 at r2500, 1.0 at r500 (logged)",
                             "class scatters SPARC .20 / ETG .22 / S09 .12 / V06 .08 (empirical)"],
               field_ETG_rule="rho_env < 0.5 mechanical (residual group centrals noted: N3607, N1407)"),
          open(os.path.join(DERIVED, "t3_3_straddle.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.0, 4.0))
colors = {"SPARC": ("#1f77b4", 7, .3), "ETG_field": ("#2ca02c", 26, .9), "Sun09": ("#7f4fc9", 26, .85), "V06": ("crimson", 40, .95)}
for c, (col, s, al) in colors.items():
    m = cls == c
    ax.scatter(sig[m], D[m], s=s, color=col, alpha=al, label=f"{c} (N={m.sum()})", zorder=3 if c != "SPARC" else 2)
xs = np.geomspace(20, 2000, 300)
ax.plot(xs, model(np.log(xs), *p0), "k-", lw=1.6, zorder=4, label="logistic fit")
ax.axvspan(sb_lo, sb_hi, color="k", alpha=.10, lw=0)
ax.axvline(sb_med, color="k", lw=1)
ax.text(sb_med*1.04, -.36, f"$\\sigma_b={sb_med:.0f}^{{+{sb_hi-sb_med:.0f}}}_{{-{sb_med-sb_lo:.0f}}}$ km/s\n$\\Rightarrow m = {m_med:.1f}$ eV", fontsize=8)
ax.axvspan(476, 756, color="orange", alpha=.10, lw=0)
ax.axvline(600, color="orange", ls=":", lw=1.2)
ax.text(610, .62, "frozen fiducial\n(m=8.45 eV, $\\rho$ band)", fontsize=7, color="darkorange")
ax.axhline(0, color="k", lw=.6)
ax.set(xscale="log", xlim=(20, 2000), ylim=(-.45, .78),
       xlabel=r"host-system $\sigma$ [km s$^{-1}$] (frozen proxies)", ylabel=r"$D$ [dex]")
ax.legend(loc="upper left", fontsize=7)
ax.set_title("T3.3: the straddle scan — break fit and carrier-mass inversion", fontsize=9.5)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_t3_dsigma.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/t3_3_straddle.json; updated figures/fig_t3_dsigma.png")