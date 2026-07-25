#!/usr/bin/env python
"""
T3.1 — anchors: the two plateaus of D(sigma).
Condensed side: SPARC (pre-freeze staged, anchor-only per plan) with locked a0.
Normal side + straddle seeds: Vikhlinin+06 (13 relaxed Chandra clusters/groups),
parsed at runtime from the staged e-print TeX (data/T3/vikhlinin06/mprof.tex):
  g_obs from their NFW fits (M2500 normalization, r_s = r500/c500),
  g_bar from their Table-2 gas-density models (rho_g = 1.624 m_p sqrt(np*ne)) + 1e12 Msun BCG.
D per system over the frozen g_bar window; sigma from the frozen proxies
(sigma = v_flat/sqrt(2) for discs; 400.4 sqrt(T_spec/keV) for clusters).
"""
import os, sys, re, json, zipfile, io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import enforce, ROOT, DERIVED

T3 = os.path.join(ROOT, "data/T3")
STAGED = [os.path.join(T3, "vikhlinin06/mprof.tex"),
          os.path.join(T3, "accept/table1.dat"), os.path.join(T3, "accept/table5.dat")]
freezes, _ = enforce(staged=STAGED, allow_dirty_paths=("t3_1_anchors.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
A0_EMP = FZ["statistic"]["g_RAR_fork_empirical"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
print(f"T3 guard OK; frozen a0={A0}, window=[{GLO},{GHI}], H0={H0}")

G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19; MP = 1.6726e-27
def g_rar(gbar, a0): return gbar / (1 - np.exp(-np.sqrt(gbar / a0)))

def D_of(r_m, gobs, gbar, a0=A0):
    m = (gbar >= GLO) & (gbar <= GHI) & (gobs > 0)
    if m.sum() < 3: return None
    return float(np.mean(np.log10(gobs[m]) - np.log10(g_rar(gbar[m], a0))))

# ---------- Vikhlinin+06 ----------
tex = open(STAGED[0]).read().replace(r"\pz", " ").replace("~", " ")
blocks = re.findall(r"\\startdata(.*?)\\enddata", tex, re.S)
def rows(block):
    out = []
    for line in re.split(r"\\\\", block):
        line = line.replace("\\dotfill", "").replace("$", "")
        line = re.sub(r"\\mcc\{[^}]*\}", "nan", line).replace(r"\nodata", "nan")
        line = re.sub(r"\\,", "", line)
        cells = [c.strip() for c in line.split("&")]
        if len(cells) < 3: continue
        name = re.sub(r"[\s{}\\]", "", cells[0])
        vals = []
        for c in cells[1:]:
            mm = re.match(r"^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", c.replace("+", "", 1) if c.strip().startswith("+") and len(c.strip()) == 1 else c)
            vals.append(float(mm.group(1)) if mm else np.nan)
        out.append((name, vals))
    return out

t_sample = {n: v for n, v in rows(blocks[0])}          # z, rmin, rdet
t_dens   = {n: v for n, v in rows(blocks[1])}          # rdet n0 rc rs a b eps n02 rc2 b2
t_mass   = {n: v for n, v in rows(blocks[2])}          # r500 Tspec Tmg c500 M2500 M500 fg...
rho_c0 = 3 * (H0 * 1e3 / (KPC * 1e3)) ** 2 / (8 * np.pi * G)   # kg/m^3

def npne(r_kpc, p):
    n0, rc, rs, al, be, ep, n02, rc2, b2 = p
    t1 = (n0 * 1e-3) ** 2 * (r_kpc / rc) ** (-al) / (1 + r_kpc**2 / rc**2) ** (3*be - al/2) \
         / (1 + r_kpc**3 / rs**3) ** (ep / 3)
    t2 = 0.0 if np.isnan(n02) else (n02 * 1e-1) ** 2 / (1 + r_kpc**2 / rc2**2) ** (3*b2)
    return t1 + t2   # cm^-6

V, skipped = [], []
for name, sv in t_sample.items():
    if name not in t_dens or name not in t_mass: continue
    z, rmin, rdet = sv[0], sv[1], sv[2]
    dm = t_dens[name][1:10]
    r500, Tspec, Tmg, c500, M2500 = (t_mass[name] + [np.nan]*9)[0:5]
    if any(np.isnan(x) for x in (r500, c500, M2500, Tspec)):
        skipped.append(name); continue
    rho_cz = rho_c0 * (OM * (1+z)**3 + 1 - OM)
    r2500 = (3 * M2500 * 1e14 * MSUN / (4*np.pi * 2500 * rho_cz)) ** (1/3) / KPC   # kpc
    rs_nfw = r500 / c500
    mu = lambda x: np.log(1+x) - x/(1+x)
    Mnfw = lambda r: M2500 * 1e14 * MSUN * mu(r/rs_nfw) / mu(r2500/rs_nfw)
    rg = np.geomspace(max(rmin, 8), rdet, 240)
    rho_g = 1.624 * MP * np.sqrt(npne(rg, dm)) * 1e3   # kg/m^3 (cm^-3 -> m^-3 factor 1e6 inside sqrt of cm^-6 => 1e3... )
    rho_g = 1.624 * MP * np.sqrt(npne(rg, dm) * 1e12)  # np*ne in cm^-6 -> m^-6
    Mg = np.concatenate([[0], np.cumsum(0.5*(rho_g[1:]*4*np.pi*(rg[1:]*KPC)**2 + rho_g[:-1]*4*np.pi*(rg[:-1]*KPC)**2) * np.diff(rg)*KPC)])
    Mbar = Mg + 1e12 * MSUN
    gobs = G * Mnfw(rg) / (rg*KPC)**2
    gbar = G * Mbar / (rg*KPC)**2
    D  = D_of(rg, gobs, gbar); De = D_of(rg, gobs, gbar, A0_EMP)
    sig = 400.4 * np.sqrt(Tspec)
    if D is not None:
        V.append(dict(name=name, z=z, sigma=round(sig,0), Tspec=Tspec, D=round(D,3), D_empfork=round(De,3),
                      r500=r500, M2500e14=M2500))
print(f"V06 systems used: {len(V)}; skipped (incomplete NFW): {skipped}")

# ---------- SPARC ----------
zf = zipfile.ZipFile(os.path.join(ROOT, "data/T2/sparc/Rotmod_LTG.zip"))
S = []
for nm in zf.namelist():
    if not nm.endswith("_rotmod.dat"): continue
    try:
        arr = np.loadtxt(io.BytesIO(zf.read(nm)), comments="#")
    except Exception: continue
    if arr.ndim != 2 or arr.shape[0] < 6 or arr.shape[1] < 6: continue
    r, vobs, _, vgas, vdisk, vbul = (arr[:, i] for i in range(6))
    r_m = r * KPC; km = 1e3
    gobs = (vobs*km)**2 / r_m
    gbar = (vgas*km*np.abs(vgas*km) + 0.5*(vdisk*km)**2 + 0.7*(vbul*km)**2) / r_m
    ok = gbar > 0
    D = D_of(r_m[ok], gobs[ok], gbar[ok]); De = D_of(r_m[ok], gobs[ok], gbar[ok], A0_EMP)
    vflat = float(np.mean(vobs[-3:]))
    if D is not None and vflat > 0:
        S.append(dict(name=nm.split("_rotmod")[0], sigma=round(vflat/np.sqrt(2),1), D=round(D,3), D_empfork=round(De,3)))

sD  = np.array([x["D"] for x in S]); ssig = np.array([x["sigma"] for x in S])
p1m = ssig < 200
P1_mean, P1_scatter = float(np.mean(sD[p1m])), float(np.std(sD[p1m]))
vD  = np.array([x["D"] for x in V]); vsig = np.array([x["sigma"] for x in V])
hi  = vsig > 1000
P2_mean = float(np.mean(vD[hi])) if hi.any() else None

print("\n=== T3.1 ANCHORS ===")
print(f"P1 condensed floor (SPARC, sigma<200, N={p1m.sum()}): mean D = {P1_mean:+.3f} dex "
      f"(scatter {P1_scatter:.3f}) — frozen threshold |D|<0.10: {'PASS' if abs(P1_mean)<0.10 else 'FAIL'}")
if P2_mean is not None:
    print(f"P2 normal excess (V06, sigma>1000, N={int(hi.sum())}): mean D = {P2_mean:+.3f} dex "
          f"— frozen threshold D>0.15: {'PASS' if P2_mean>0.15 else 'FAIL'}")
print("\nper-cluster (V06):")
for x in sorted(V, key=lambda y: y["sigma"]):
    print(f"  {x['name']:14s} T={x['Tspec']:5.2f} keV  sigma={x['sigma']:6.0f}  D={x['D']:+.3f}  (emp-fork {x['D_empfork']:+.3f})")

json.dump(dict(frozen_against=FZ["provenance"], P1=dict(mean=P1_mean, scatter=P1_scatter, N=int(p1m.sum())),
               P2=dict(mean=P2_mean, N=int(hi.sum())), sparc_note="pre-freeze staged, anchor-only",
               bcg_note="1e12 Msun nominal BCG added to g_bar (V06)", v06=V, sparc_n=len(S)),
          open(os.path.join(DERIVED, "t3_1_anchors.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.scatter(ssig, sD, s=8, alpha=.35, color="#1f77b4", label=f"SPARC discs (N={len(S)})")
ax.scatter(vsig, vD, s=42, color="crimson", zorder=3, label="V06 clusters/groups (NFW vs gas+BCG)")
for x in V: ax.annotate(x["name"], (x["sigma"], x["D"]), fontsize=5.5, xytext=(2,3), textcoords="offset points")
kern = FZ["kernel"]["sigma_crit_km_s"]
ax.axvspan(kern["rho_band_at_m_fid"][0], kern["rho_band_at_m_fid"][1], color="k", alpha=.08, lw=0)
ax.axvline(600, color="k", ls=":", lw=1); ax.text(608, .58, "frozen $\\sigma_{\\rm crit}$ fiducial\n(+$\\rho$ band)", fontsize=7)
ax.axvspan(kern["outer_envelope"][0], kern["outer_envelope"][1], color="orange", alpha=.05, lw=0)
ax.axhline(0, color="k", lw=.6); ax.axhline(0.15, color="crimson", lw=.6, ls="--")
ax.set(xscale="log", xlabel=r"$\sigma$ [km s$^{-1}$] (frozen proxies)", ylabel=r"$D$ [dex]",
       xlim=(15, 2000), ylim=(-.45, .75))
ax.legend(loc="upper left", fontsize=7.5)
ax.set_title("T3.1 anchors: the two plateaus (frozen statistic, locked $a_0$)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_t3_dsigma.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/t3_1_anchors.json, figures/fig_t3_dsigma.png")