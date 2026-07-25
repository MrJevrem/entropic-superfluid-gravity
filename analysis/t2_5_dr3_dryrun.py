#!/usr/bin/env python
"""
T2.5 phase 2 — DR3 dry-run (plan §7): both contested estimators on the El-Badry+21 eDR3
wide-binary catalog, under the FROZEN adjudication rules (sample cuts, e-prior forks, triple
handling, class thresholds). Purpose: pipeline validation + baseline landscape. NO branch
adjudication from eDR3 — phase 3 (DR4, Dec 2026) runs the same code on the new catalog.

Estimator A (Chae-style deprojection): data v-tilde distribution per r_p bin vs matched
Newtonian MC; delta-v-tilde = <v>_data/<v>_MC - 1.
Estimator B (Banik-style forward model): fit alpha = G_eff/G by multinomial likelihood on the
v-tilde histogram, triples nuisance f_t marginalized over the frozen grid; per e-prior fork.
Validation: high-acceleration control bin (0.2-1 kAU) must return delta~0, alpha~1.

[P] logged approximations: mass from M_G via Pecaut-Mamajek-style interpolation; triple boost
u~U(0.3,0.8) on v; pm-error correlations neglected; single-epoch projected statistic only.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from t2_guard import enforce, DERIVED, ROOT

CAT = os.path.join(ROOT, "data/T2/wb/all_columns_catalog.fits.gz")
freezes, _ = enforce(staged=[CAT])
fz = freezes[-1]
RULES = fz["wb_adjudication"]
rng = np.random.default_rng(20260724)

# ---------------- load + frozen cuts ----------------
cols = ["parallax1", "parallax2", "parallax_over_error1", "parallax_over_error2",
        "pmra1", "pmra2", "pmdec1", "pmdec2", "pmra_error1", "pmra_error2",
        "pmdec_error1", "pmdec_error2", "ruwe1", "ruwe2", "phot_g_mean_mag1",
        "phot_g_mean_mag2", "sep_AU", "binary_type", "R_chance_align"]
with fits.open(CAT, memmap=False) as h:
    d = h[1].data
    A = {c: np.asarray(d[c]) for c in cols}
    del d
n0 = len(A["sep_AU"])
plx = 0.5 * (A["parallax1"] + A["parallax2"])
dist_pc = 1000.0 / np.clip(plx, 1e-3, None)
MG1 = A["phot_g_mean_mag1"] + 5 * np.log10(np.clip(A["parallax1"], 1e-3, None)) - 10
MG2 = A["phot_g_mean_mag2"] + 5 * np.log10(np.clip(A["parallax2"], 1e-3, None)) - 10
mask = ((A["binary_type"] == "MSMS") & (A["parallax_over_error1"] > 20)
        & (A["parallax_over_error2"] > 20) & (plx > 5.0)
        & (A["ruwe1"] < 1.4) & (A["ruwe2"] < 1.4) & (A["R_chance_align"] < 0.1)
        & (A["phot_g_mean_mag1"] < 18) & (A["phot_g_mean_mag2"] < 18))
print(f"catalog {n0} pairs -> frozen-cut sample {mask.sum()}")

# mass from M_G [P]: Pecaut-Mamajek-style MS anchor points
MGgrid = np.array([1.4, 2.6, 4.0, 4.8, 5.6, 6.7, 8.0, 9.5, 11.0, 13.0, 15.0])
Mgrid = np.array([2.0, 1.5, 1.1, 1.0, 0.9, 0.75, 0.60, 0.45, 0.30, 0.15, 0.09])
def mass_of(MG): return np.interp(np.clip(MG, MGgrid[0], MGgrid[-1]), MGgrid, Mgrid)
Mtot = mass_of(MG1) + mass_of(MG2)

dmu = np.hypot(A["pmra1"] - A["pmra2"], A["pmdec1"] - A["pmdec2"])          # mas/yr
sdmu = np.sqrt(A["pmra_error1"] ** 2 + A["pmra_error2"] ** 2
               + A["pmdec_error1"] ** 2 + A["pmdec_error2"] ** 2) / np.sqrt(2)
vp = 4.74047e-3 * dmu * dist_pc                                              # km/s
svp = 4.74047e-3 * sdmu * dist_pc
vc = 29.785 * np.sqrt(Mtot / np.clip(A["sep_AU"], 1, None))                  # km/s at r_p
vt = vp / vc
good = mask & (svp < 0.15 * vc) & np.isfinite(vt)
print(f"after velocity-precision cut: {good.sum()}")

# ---------------- Newtonian MC (scale-free v-tilde distributions per e-prior) ----------------
def newton_mc(e_prior, n=400_000):
    if e_prior == "thermal":
        e = np.sqrt(rng.random(n))                       # f(e)=2e
    else:                                                # Hwang+22-like superthermal, alpha=1.2
        e = rng.random(n) ** (1 / 2.2)                   # p(e) ∝ e^1.2
    Mano = rng.random(n) * 2 * np.pi
    E = Mano.copy()
    for _ in range(12):
        E -= (E - e * np.sin(E) - Mano) / (1 - e * np.cos(E))
    cosE, sinE = np.cos(E), np.sin(E)
    b = np.sqrt(1 - e ** 2)
    r = np.stack([cosE - e, b * sinE], 1)                # a=1, GM=1
    v = np.stack([-sinE, b * cosE], 1) / (1 - e * cosE)[:, None]
    w = rng.random(n) * 2 * np.pi
    cw, sw = np.cos(w), np.sin(w)
    R = np.stack([np.stack([cw, -sw], 1), np.stack([sw, cw], 1)], 1)
    r2 = np.einsum("nij,nj->ni", R, r); v2 = np.einsum("nij,nj->ni", R, v)
    ci = rng.random(n) * 2 - 1                            # cos(i) uniform
    r3 = np.stack([r2[:, 0], r2[:, 1] * ci], 1)           # project (incline about x-axis)
    v3 = np.stack([v2[:, 0], v2[:, 1] * ci], 1)
    s = np.linalg.norm(r3, axis=1); vsky = np.linalg.norm(v3, axis=1)
    return vsky / np.sqrt(1 / s)                          # v-tilde (scale-free)

MC = {p: newton_mc(p) for p in ("thermal", "hwang")}

# ---------------- estimators per bin ----------------
BINS = {"control_0.2-1kAU": (200, 1000), "2-5kAU": (2000, 5000),
        "5-10kAU": (5000, 10000), "10-30kAU": (10000, 30000)}
HEDGES = np.linspace(0, 2.5, 26)
def chae_delta(vt_dat, mc):
    m = vt_dat[(vt_dat > 0) & (vt_dat < 2.5)]
    mcm = mc[(mc > 0) & (mc < 2.5)]
    boots = []
    for _ in range(300):
        boots.append(np.mean(rng.choice(m, len(m))) / np.mean(mcm) - 1)
    return float(np.mean(m) / np.mean(mcm) - 1), float(np.std(boots))
def banik_alpha(vt_dat, mc):
    m = vt_dat[(vt_dat > 0) & (vt_dat < 2.5)]
    obs, _ = np.histogram(m, HEDGES)
    best = None
    for alpha in np.linspace(0.7, 1.8, 56):
        for ft in (0.0, 0.05, 0.10, 0.15):
            sim = np.sqrt(alpha) * mc
            nt = int(ft * len(sim))
            if nt:
                boost = np.ones(len(sim)); idx = rng.choice(len(sim), nt, replace=False)
                boost[idx] = 1 + rng.uniform(0.3, 0.8, nt)
                sim = sim * boost
            pm, _ = np.histogram(sim[(sim > 0) & (sim < 2.5)], HEDGES)
            pm = pm / max(pm.sum(), 1) + 1e-12
            ll = float(np.sum(obs * np.log(pm)))
            if best is None or ll > best[0]: best = (ll, alpha, ft)
    return best[1], best[2]

results = {"rules": RULES, "n_sample": int(good.sum()), "bins": {}}
print(f"\n{'bin':16}{'n':>6}  {'e-prior':8}{'dv (Chae)':>16}{'alpha (Banik)':>14}{'f_t':>5}")
for bname, (lo, hi) in BINS.items():
    sel = good & (A["sep_AU"] >= lo) & (A["sep_AU"] < hi)
    vt_dat = vt[sel]
    results["bins"][bname] = {"n": int(sel.sum())}
    for prior in ("thermal", "hwang"):
        dv, dverr = chae_delta(vt_dat, MC[prior])
        al, ft = banik_alpha(vt_dat, MC[prior])
        results["bins"][bname][prior] = dict(delta_v=round(dv, 4), delta_v_err=round(dverr, 4),
                                             alpha=round(al, 3), f_t=ft)
        print(f"{bname:16}{sel.sum():6d}  {prior:8}{dv:+9.3f}±{dverr:.3f}"
              f"{al:14.3f}{ft:5.2f}")

# ---------------- validation + landscape ----------------
ctrl = results["bins"]["control_0.2-1kAU"]
val_ok = all(abs(ctrl[p]["delta_v"]) < 0.05 and abs(ctrl[p]["alpha"] - 1) < 0.1
             for p in ("thermal", "hwang"))
results["validation"] = dict(
    control_bin_newtonian=bool(val_ok),
    note="control bin (g >> a0) must give dv~0, alpha~1 under both e-priors")
results["dry_run_status"] = (
    "PIPELINE " + ("VALIDATED" if val_ok else "NOT VALIDATED") + " on eDR3. Published-result "
    "reproduction (Chae 2023 boost; Banik 2024 alpha~1) is qualitative only at this stage — "
    "full per-paper reproduction remains open engineering before Dec 2026. NO branch "
    "adjudication from this dry-run (frozen: phase 3 = DR4 with identical code).")
json.dump(results, open(os.path.join(DERIVED, "t2_5_results.json"), "w"), indent=1)
print("\nvalidation:", results["validation"])
print("saved derived/t2_5_results.json")
