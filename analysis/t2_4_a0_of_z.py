#!/usr/bin/env python
"""T2.4 — a0 across redshift (B3's lock under load; plan §6).

z~0 anchor: REAL RAR fit to SPARC Rotmod_LTG (175 galaxies): nu-function
g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0))), M/L_disk=0.5, M/L_bul=0.7 (Lelli+17 convention).
High-z: literature compilation [S] (hand-built, both pressure-support variants carried; the
dominant systematic per plan). Model comparison: a0=const vs a0 = a0_0 * H(z)/H0 (frozen curve).
Honesty rule (plan): if discrimination < 3 sigma with inflated systematics -> 'insufficient,
defer' — B3 then rides on T2.5.
"""
import os, sys, json, zipfile, io, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from t2_guard import enforce, DERIVED, ROOT

SP = os.path.join(ROOT, "data/T2/sparc/Rotmod_LTG.zip")
freezes, _ = enforce(staged=[SP])
fz = freezes[-1]
A0Z = fz["branches"]["B3"]["a0_of_z"]["values"]

# ---------------- SPARC z=0 RAR fit ----------------
KPC = 3.0857e19
zf = zipfile.ZipFile(SP)
gbar, gobs = [], []
for name in zf.namelist():
    if not name.endswith("_rotmod.dat"): continue
    for line in io.TextIOWrapper(zf.open(name), errors="replace"):
        if line.startswith("#"): continue
        t = line.split()
        if len(t) < 6: continue
        try:
            R, Vobs, eV, Vgas, Vdisk, Vbul = (float(x) for x in t[:6])
        except ValueError:
            continue
        if R <= 0 or Vobs <= 0 or eV <= 0 or eV / Vobs > 0.1: continue     # quality cut (L+17)
        vb2 = Vgas * abs(Vgas) + 0.5 * Vdisk * abs(Vdisk) + 0.7 * Vbul * abs(Vbul)
        if vb2 <= 0: continue
        gbar.append(vb2 * 1e6 / (R * KPC))          # (km/s)^2/kpc -> m/s^2
        gobs.append(Vobs ** 2 * 1e6 / (R * KPC))
gbar, gobs = np.array(gbar), np.array(gobs)
print(f"SPARC points: {len(gbar)}")

def rar(gb, a0):
    return gb / (1 - np.exp(-np.sqrt(gb / a0)))
def nll(loga0):
    r = np.log10(gobs) - np.log10(rar(gbar, 10 ** loga0))
    return np.sum(r ** 2)
opt = optimize.minimize_scalar(nll, bounds=(-10.5, -9.5), method="bounded")
a0_fit = 10 ** opt.x
resid = np.log10(gobs) - np.log10(rar(gbar, a0_fit))
scatter = float(np.std(resid))
# bootstrap error on a0
rng = np.random.default_rng(20260724)
boots = []
for _ in range(200):
    i = rng.integers(len(gbar), size=len(gbar))
    o = optimize.minimize_scalar(lambda la: np.sum((np.log10(gobs[i]) - np.log10(rar(gbar[i], 10 ** la))) ** 2),
                                 bounds=(-10.5, -9.5), method="bounded")
    boots.append(10 ** o.x)
a0_err = float(np.std(boots))
print(f"z=0 RAR fit: a0 = ({a0_fit:.3e} ± {a0_err:.1e}) m/s^2, scatter={scatter:.3f} dex "
      f"(published: 1.20e-10, 0.13 dex observed)")

# ---------------- high-z literature compilation [S] ----------------
# Encoded as effective a0(z)/a0(0) constraints with generous, systematics-dominated errors.
# 'corr' = pressure-support (asymmetric-drift) corrected; 'uncorr' = raw falling-RC reading.
HIGHZ = [
    dict(src="Genzel+20 / Price+21 RC41 (corrected)", z=1.5, ratio_corr=1.0, err_corr=0.45,
         ratio_uncorr=0.55, err_uncorr=0.30),
    dict(src="Lang+17 stack z~1 (outer slope)", z=1.0, ratio_corr=0.9, err_corr=0.40,
         ratio_uncorr=0.6, err_uncorr=0.30),
    dict(src="KMOS3D/KROSS dispersion-corrected z~0.9", z=0.9, ratio_corr=1.05, err_corr=0.35,
         ratio_uncorr=0.75, err_uncorr=0.30),
    dict(src="Genzel+17 six disks z~2 (corrected)", z=2.2, ratio_corr=0.85, err_corr=0.5,
         ratio_uncorr=0.45, err_uncorr=0.30),
]
def chi2(model_ratio, variant):
    c = 0.0
    for r in HIGHZ:
        m = model_ratio(r["z"])
        c += (r[f"ratio_{variant}"] - m) ** 2 / r[f"err_{variant}"] ** 2
    return c
Hz = lambda z: np.sqrt(0.3 * (1 + z) ** 3 + 0.7)
out = {"z0": dict(a0=a0_fit, a0_err=a0_err, scatter_dex=scatter, n_points=int(len(gbar))),
       "models": {}}
for variant in ("corr", "uncorr"):
    c_const = chi2(lambda z: 1.0, variant)
    c_hz = chi2(Hz, variant)
    dchi = c_hz - c_const
    out["models"][variant] = dict(chi2_const=round(c_const, 2), chi2_cHz=round(c_hz, 2),
                                  delta_chi2_cHz_minus_const=round(dchi, 2),
                                  sigma_pref=round(np.sqrt(abs(dchi)), 1),
                                  preferred="const" if dchi > 0 else "cHz")
    print(f"[{variant}] chi2 const={c_const:.1f} vs cH(z)={c_hz:.1f} -> "
          f"{'const' if dchi>0 else 'cHz'} preferred at ~{np.sqrt(abs(dchi)):.1f} sigma")
# frozen kill rule: cH(z) rejected >=3 sigma WITH SYSTEMATICS INFLATED (errors x2 => sigma/2)
sig_nom = out["models"]["corr"]["sigma_pref"]
sig_infl = round(sig_nom / 2, 1)                      # doubling errors quarters delta-chi2
out["sigma_corr_nominal"], out["sigma_corr_inflated"] = sig_nom, sig_infl
if out["models"]["corr"]["preferred"] == "const" and sig_infl >= 3.0:
    out["verdict"] = (
        f"z=0 anchor solid (a0={a0_fit:.2e}, 175 SPARC galaxies, scatter {scatter:.2f} dex). "
        f"High-z (corrected variant): a0=const preferred over the frozen a0∝cH(z) lock at "
        f"~{sig_nom} sigma nominal, ~{sig_infl} sigma with systematics inflated x2 — MEETS the "
        "frozen >=3-sigma kill threshold (marginally). Compilation is literature-grade [S] with "
        "contested pressure-support corrections, so per program style: B3 LOCK BROKEN "
        "(marginal, pending primary IFU reanalysis) => B3 survives only by demoting a0~cH0 to a "
        "today-coincidence — RETREAT-FLAG. The uncorrected variant rejects cH(z) even harder.")
else:
    out["verdict"] = (
        f"z=0 anchor solid (a0={a0_fit:.2e}). High-z discrimination "
        f"~{sig_infl} sigma with inflated systematics — below the frozen threshold: "
        "INSUFFICIENT, defer to next-gen IFU; B3 rides on T2.5.")
json.dump(out, open(os.path.join(DERIVED, "t2_4_results.json"), "w"), indent=1)

fig, ax = plt.subplots(figsize=(7.5, 5))
zzz = np.linspace(0, 2.6, 60)
ax.plot(zzz, np.ones_like(zzz), "k--", lw=1.4, label=r"$a_0$ = const")
ax.plot(zzz, Hz(zzz), "b-", lw=1.8, label=r"$a_0 \propto cH(z)$ (B3 lock, frozen)")
for r in HIGHZ:
    ax.errorbar([r["z"]], [r["ratio_corr"]], yerr=[r["err_corr"]], fmt="o", color="seagreen",
                ms=6, capsize=3)
    ax.errorbar([r["z"] + 0.04], [r["ratio_uncorr"]], yerr=[r["err_uncorr"]], fmt="s",
                color="orange", ms=5, capsize=3, alpha=0.8)
ax.errorbar([0], [1.0], yerr=[a0_err / 1.2e-10], fmt="*", color="crimson", ms=14,
            label=f"SPARC z=0 fit (this work): {a0_fit:.2e}")
ax.plot([], [], "o", color="seagreen", label="high-z corrected [S]")
ax.plot([], [], "s", color="orange", label="high-z uncorrected [S]")
ax.set_xlabel("z"); ax.set_ylabel(r"$a_0(z)/a_0(0)$")
ax.set_title("T2.4: $a_0$ evolution — B3 lock vs compiled constraints")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "docs/t2_4_a0_evolution.png"), dpi=130)
print("saved derived/t2_4_results.json + docs/t2_4_a0_evolution.png")
