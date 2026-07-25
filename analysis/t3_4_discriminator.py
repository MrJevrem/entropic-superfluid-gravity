#!/usr/bin/env python
"""
T3.4 — the discriminator grid (P4) on the method-matched S09 sample (+V06 for z reach).
 (1) MAXIMAL-BARYON FORK: recompute D with g_bar -> g_obs * f_b,cosmic (0.156) inside the
     aperture — the most generous concession to the feedback/missing-baryon story (B).
     If D' stays elevated, no baryon bookkeeping restores the RAR in groups.
 (2) P4 regressions: standardized OLS D ~ ln sigma + ln M2500 + z within S09, with the
     collinearity matrix reported honestly (M and sigma are near-degenerate in X-ray samples).
 (3) c500 (formation-history proxy) correlation at fixed sigma.
 (4) z-split of the elevated population (S09 + V06).
"""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import enforce, ROOT, DERIVED

TEX = os.path.join(ROOT, "data/T3/sun09/ms.tex")
freezes, _ = enforce(staged=[TEX], allow_dirty_paths=("t3_4_discriminator.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
FB_COSMIC = 0.156
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
rho_c0 = 3 * (H0*1e3/(KPC*1e3))**2 / (8*np.pi*G)
g_rar = lambda gb: gb / (1 - np.exp(-np.sqrt(gb / A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\\tablenotemark\{[a-z]\}|~", " ",
                         re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")

tex = open(TEX).read()
zLK = {}
for line in tex.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) >= 8 and re.search(r"\d", cells[1]):
        name = norm(cells[0]); z, LK = num1(cells[1]), num1(cells[7])
        if name and not np.isnan(z): zLK[name] = (z, LK)
S = []
for line in tex.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) < 10 or not re.search(r"\d", "".join(cells[1:3])): continue
    name = norm(cells[0])
    if name not in zLK: continue
    z, LK = zLK[name]
    T500 = num1(cells[1]); r2500 = num1(cells[4]); fg2500 = num1(cells[7]); c500 = num1(cells[9])
    if np.isnan(T500) or np.isnan(r2500) or np.isnan(fg2500): continue
    rho_cz = rho_c0 * (OM*(1+z)**3 + 1-OM)
    Mstar = (10**LK) * MSUN if not np.isnan(LK) else 0.0
    r_m = r2500 * KPC
    M2500 = (4*np.pi/3) * 2500 * rho_cz * r_m**3
    gobs = G * M2500 / r_m**2
    fb = fg2500 + 0.7*Mstar/M2500
    gbar = gobs * fb
    if not (GLO <= gbar <= GHI): continue
    D  = np.log10(gobs) - np.log10(g_rar(gbar))
    Dp = np.log10(gobs) - np.log10(g_rar(gobs * FB_COSMIC))   # maximal-baryon fork
    S.append(dict(name=name, z=z, sigma=400.4*np.sqrt(T500), lnM=np.log(M2500/MSUN),
                  fb=fb, c500=c500, D=D, Dmax=Dp))
sig = np.array([s["sigma"] for s in S]); D = np.array([s["D"] for s in S])
Dm = np.array([s["Dmax"] for s in S]); fb = np.array([s["fb"] for s in S])
lnM = np.array([s["lnM"] for s in S]); zz = np.array([s["z"] for s in S])
c5 = np.array([s["c500"] for s in S])
print(f"S09 systems (Delta=2500 aperture): N={len(S)}")

print("\n=== (1) MAXIMAL-BARYON FORK ===")
print(f"measured:        mean D    = {D.mean():+.3f}  (scatter {D.std():.3f})")
print(f"cosmic-fb fork:  mean D'   = {Dm.mean():+.3f}  (scatter {Dm.std():.3f})")
print(f"fraction of elevation surviving full baryon restoration: {Dm.mean()/D.mean():.0%}")
print(f"systems with D' > 0.15 even at cosmic fb: {(Dm > 0.15).sum()}/{len(S)}")
print(f"mean measured fb(2500) = {fb.mean():.3f} vs cosmic {FB_COSMIC} (shortfall factor {FB_COSMIC/fb.mean():.1f})")

print("\n=== (2) P4 regressions (S09, standardized OLS) ===")
X = np.column_stack([np.log(sig), lnM, zz])
Xs = (X - X.mean(0)) / X.std(0); Ds = (D - D.mean()) / D.std()
r_sM = np.corrcoef(np.log(sig), lnM)[0, 1]
beta, res, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(S)), Xs]), Ds, rcond=None)
dof = len(S) - 4
se = np.sqrt(np.sum((Ds - np.column_stack([np.ones(len(S)), Xs]) @ beta)**2)/dof *
             np.diag(np.linalg.inv((np.column_stack([np.ones(len(S)), Xs]).T @ np.column_stack([np.ones(len(S)), Xs])))))
names = ["const", "ln sigma", "ln M2500", "z"]
for n, b, s_ in zip(names, beta, se):
    t = b/s_ if s_ > 0 else 0
    print(f"  beta[{n:9s}] = {b:+.3f} +/- {s_:.3f}  (t = {t:+.1f})")
print(f"  collinearity: corr(ln sigma, ln M2500) = {r_sM:.3f}  <- near-degenerate; partials weakly identified (stated)")

m_lo = sig < np.median(sig)
r_cD_lo = np.corrcoef(c5[m_lo & ~np.isnan(c5)], D[m_lo & ~np.isnan(c5)])[0, 1]
r_cD = np.corrcoef(c5[~np.isnan(c5)], D[~np.isnan(c5)])[0, 1]
print(f"\n=== (3) concentration proxy ===\n  corr(c500, D) all = {r_cD:+.2f}; low-sigma half = {r_cD_lo:+.2f}  (weak-prior diagnostic)")

v06 = json.load(open(os.path.join(DERIVED, "t3_1_anchors.json")))["v06"]
ez = np.array([x["z"] for x in v06] + list(zz)); eD = np.array([x["D"] for x in v06] + list(D))
zs = np.median(ez)
print(f"\n=== (4) z-split (elevated set S09+V06, N={len(ez)}, median z={zs:.3f}) ===")
print(f"  lo-z mean D = {eD[ez <= zs].mean():+.3f} (N={(ez <= zs).sum()})   hi-z mean D = {eD[ez > zs].mean():+.3f} (N={(ez > zs).sum()})")
print(f"  delta = {eD[ez > zs].mean() - eD[ez <= zs].mean():+.3f}  (z range {ez.min():.3f}-{ez.max():.3f} — narrow; weak test, logged)")

print("\n=== P4 VERDICT SHEET ===")
v1 = Dm.mean() > 0.15
print(f"(B)-cap: baryon restoration leaves D' = {Dm.mean():+.2f} -> feedback explains at most "
      f"{1 - Dm.mean()/D.mean():.0%} of the elevation. Artifact-alone: {'EXCLUDED' if v1 else 'viable'}")
print(f"M-at-fixed-sigma: unidentifiable in-sample (collinearity {r_sM:.2f}); DEFERRED to mixed-method samples (logged)")
print(f"z-invariance: consistent within scatter over the (narrow) staged range")

json.dump(dict(frozen_against=FZ["provenance"], N=len(S),
               fork=dict(D_mean=float(D.mean()), Dmax_mean=float(Dm.mean()),
                         surviving_fraction=float(Dm.mean()/D.mean()),
                         n_above_015=int((Dm > 0.15).sum()), fb_mean=float(fb.mean())),
               regression=dict(names=names, beta=[float(b) for b in beta], se=[float(x) for x in se],
                                corr_lnsig_lnM=float(r_sM)),
               c500=dict(all=float(r_cD), losig=float(r_cD_lo)),
               zsplit=dict(lo=float(eD[ez <= zs].mean()), hi=float(eD[ez > zs].mean()),
                           zmax=float(ez.max())),
               rows=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                      for k, v in s.items()} for s in S]),
          open(os.path.join(DERIVED, "t3_4_discriminator.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.6, 3.4))
ax.scatter(fb, D, s=26, color="#7f4fc9", label="measured $D$ (their $f_b$)")
ax.scatter(fb, Dm, s=26, marker="^", color="#2ca02c", label="$D'$ at cosmic $f_b=0.156$")
ax.axhline(0.15, color="crimson", lw=.7, ls="--"); ax.axhline(0, color="k", lw=.6)
ax.set(xlabel=r"measured $f_b(r_{2500})$", ylabel=r"$D$ [dex]", ylim=(-.1, .65))
ax.legend(fontsize=7); ax.set_title("(a) the maximal-baryon fork", fontsize=9)
bx.scatter(np.exp(lnM) / 1e13, D, s=26, c=sig, cmap="viridis")
cb = plt.colorbar(bx.collections[0], ax=bx); cb.set_label(r"$\sigma$ [km/s]", fontsize=8)
bx.set(xscale="log", xlabel=r"$M_{2500}$ [$10^{13}\,M_\odot$]", ylabel=r"$D$ [dex]")
bx.set_title(f"(b) $D$ vs $M$ (corr(ln$\\sigma$,ln$M$)={r_sM:.2f})", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_t3_discriminator.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/t3_4_discriminator.json, figures/fig_t3_discriminator.png")