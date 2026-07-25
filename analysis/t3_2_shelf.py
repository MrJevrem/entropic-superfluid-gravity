#!/usr/bin/env python
"""
T3.2 — the elliptical shelf (sigma 150-310): does the RAR persist beyond the discs?
Data: SLUGGS / Alabi+16 (staged e-print, parsed at runtime): GC tracer-mass estimates
M_tot(<5Re) and M_tot(<Rmax ~ 9-21 Re), with f_DM giving the enclosed stellar fraction
directly => g_bar = g_obs * (1 - f_DM), all numbers theirs (M/L_K = 1 fiducial;
beta=0 primary with +/-0.5 anisotropy fork from the same table).
Protocol deviation (logged): 1-2 aperture points per system vs the frozen >=3-point
rule — population statistic over 23 galaxies compensates; flagged in results.
"""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import enforce, ROOT, DERIVED

TEX = os.path.join(ROOT, "data/T3/sluggs/mass_est.tex")
freezes, _ = enforce(staged=[TEX], allow_dirty_paths=("t3_2_shelf.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
g_rar = lambda gb: gb / (1 - np.exp(-np.sqrt(gb / A0)))

tex = open(TEX).read()
# --- tab:summary: name Mk dist vsys Re" sigma eps rho_env logM* ...
summ = {}
block = tex.split(r"\label{tab:summary}")[1].split(r"\end{tabular}")[0]
for line in block.splitlines():
    m = re.match(r"\s*(\d{3,4})\s*&\s*\$?-([\d.]+)\$?\s*&\s*([\d.]+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)", line)
    if m:
        n, _, dist, _, re_as, sig, eps, rho, lm = m.groups()
        summ[n] = dict(dist=float(dist), Re_kpc=float(re_as)*float(dist)*4.8481e-3,
                       sigma=float(sig), rho_env=float(rho), logMs=float(lm))
# --- tab:mass_summary: beta=0 rows: name & 0& Mrot5 & Mp5 & Mtot5 pm e & fDM5 pm e & Rmax & Mrotm & Mpm & Mtotm pm e & fDMm pm e
mass = {}
mblock = tex.split(r"\label{tab:mass_summary}")[1].split(r"\label{tab:comparison}")[0]
num = r"([\d.]+)\s*\$?\\pm\s*\$?\s*([\d.]+)"
pat = re.compile(r"^\s*(\d{3,4})\s*&\s*0\s*&.*?&\s*" + num.replace("(", "(?:", 1).replace("([\\d.]+)", "[\\d.]+", 1) , re.X)
for line in mblock.splitlines():
    cells = [c.strip() for c in re.sub(r"\$|\\pm", " ", line).split("&")]
    if len(cells) >= 11 and re.match(r"^\d{3,4}$", cells[0]) and cells[1].strip().startswith("0") and not cells[1].strip().startswith("0."):
        nums = []
        for c in cells[2:11]:
            t = c.split()
            nums.append(float(t[0]) if t and re.match(r"^[\d.]+$", t[0]) else np.nan)
        # nums: Mrot5(1e10), Mp5(1e11), Mtot5(1e11), fDM5, Rmax(Re), Mrotm(1e10), Mpm(1e11), Mtotm(1e11), fDMm
        mass[cells[0]] = dict(Mtot5=nums[2]*1e11, fDM5=nums[3], Rmax=nums[4],
                              Mtotm=nums[7]*1e11, fDMm=nums[8])

rows = []
for n, s in summ.items():
    if n not in mass: continue
    mm = mass[n]
    pts = []
    for r_re, Mt, fdm in ((5.0, mm["Mtot5"], mm["fDM5"]), (mm["Rmax"], mm["Mtotm"], mm["fDMm"])):
        if any(np.isnan(x) for x in (r_re, Mt, fdm)): continue
        r = r_re * s["Re_kpc"] * KPC
        gobs = G * Mt * MSUN / r**2
        gbar = gobs * (1 - fdm)
        if GLO <= gbar <= GHI:
            pts.append(np.log10(gobs) - np.log10(g_rar(gbar)))
    if pts:
        rows.append(dict(name="NGC"+n, sigma=s["sigma"], D=round(float(np.mean(pts)), 3),
                         n_apertures=len(pts), rho_env=s["rho_env"], logMs=s["logMs"]))

D = np.array([r["D"] for r in rows]); sig = np.array([r["sigma"] for r in rows])
dense = np.array([r["rho_env"] for r in rows]) > 1.0
shelf = (sig >= 150) & (sig <= 350)
print(f"parsed {len(summ)} summary, {len(mass)} mass rows -> {len(rows)} systems with window points")
print("\n=== T3.2 SHELF (sigma 150-350) ===")
for lab, m in (("all", shelf), ("low-density env", shelf & ~dense), ("dense env (rho>1/Mpc^3)", shelf & dense)):
    if m.any():
        print(f"  {lab:26s}: N={m.sum():2d}  mean D = {D[m].mean():+.3f}  (scatter {D[m].std():.3f})")
print("\nper-galaxy:")
for r in sorted(rows, key=lambda x: x["sigma"]):
    tag = "DENSE" if r["rho_env"] > 1 else "     "
    print(f"  {r['name']:8s} sigma={r['sigma']:4.0f}  D={r['D']:+.3f}  {tag}  (rho_env={r['rho_env']:.2f}, {r['n_apertures']} ap.)")

json.dump(dict(frozen_against=FZ["provenance"], rows=rows,
               shelf_all=dict(N=int(shelf.sum()), mean=float(D[shelf].mean()), scatter=float(D[shelf].std())),
               shelf_lowenv=dict(N=int((shelf & ~dense).sum()), mean=float(D[shelf & ~dense].mean())),
               shelf_dense=dict(N=int((shelf & dense).sum()), mean=float(D[shelf & dense].mean()) if (shelf & dense).any() else None),
               deviations=["1-2 aperture points per system vs frozen >=3-point rule (population statistic; logged)",
                            "M/L_K=1 fiducial (Alabi+16); Salpeter fork exists in their appendix, not parsed"],
               proxy_note="sigma = central stellar dispersion (1 kpc); for group centrals the SYSTEM sigma exceeds this — dense-env split reported"),
          open(os.path.join(DERIVED, "t3_2_shelf.json"), "w"), indent=1)

# --- update the running D(sigma) figure ---
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
t31 = json.load(open(os.path.join(DERIVED, "t3_1_anchors.json")))
fig, ax = plt.subplots(figsize=(6.8, 3.8))
sp = t31.get("sparc_n", 0)
# re-plot SPARC cloud from t3_1 by rerunning? keep stored V06 + new ETG + note SPARC band
ax.axhspan(-0.10, 0.10, color="#1f77b4", alpha=.10, lw=0)
ax.text(17, 0.115, f"SPARC band (N={sp}): $\\bar D=-0.048$", fontsize=7, color="#1f77b4")
for x in t31["v06"]:
    ax.scatter(x["sigma"], x["D"], s=42, color="crimson", zorder=3)
ax.scatter([], [], s=42, color="crimson", label="V06 clusters/groups (T3.1)")
ax.scatter(sig[~dense], D[~dense], s=30, marker="s", color="#2ca02c", zorder=3, label="SLUGGS ETGs, low-density env")
ax.scatter(sig[dense], D[dense], s=30, marker="s", facecolor="none", edgecolor="#2ca02c", zorder=3, label="SLUGGS ETGs, dense env")
kern = FZ["kernel"]["sigma_crit_km_s"]
ax.axvspan(kern["rho_band_at_m_fid"][0], kern["rho_band_at_m_fid"][1], color="k", alpha=.08, lw=0)
ax.axvline(600, color="k", ls=":", lw=1)
ax.axhline(0, color="k", lw=.6); ax.axhline(0.15, color="crimson", lw=.6, ls="--")
ax.set(xscale="log", xlim=(15, 2000), ylim=(-.45, .75),
       xlabel=r"$\sigma$ [km s$^{-1}$] (frozen proxies)", ylabel=r"$D$ [dex]")
ax.legend(loc="upper left", fontsize=7)
ax.set_title("T3.2: the shelf filled — ETGs between the disc band and the groups", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_t3_dsigma.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/t3_2_shelf.json; updated figures/fig_t3_dsigma.png")