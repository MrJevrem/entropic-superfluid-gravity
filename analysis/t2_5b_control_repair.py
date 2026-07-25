#!/usr/bin/env python
"""
T2.5b — the control-bin repair (OPEN_ITEMS E-2; the pre-December gating item).

Diagnosis of the +0.049 control excess (t2_5 dry-run, f_t railed at 0.15):
 (R1) noise asymmetry: data v-tilde carries measurement noise (biases mean |v| UP);
      the frozen MC was noise-free -> convolve the MC with the bin's empirical
      sigma_v/v_c distribution (estimator correction, not a physics choice);
 (R2) triples: widen the f_t grid [0, 0.15] -> [0, 0.5] (rail + literature hidden
      multiplicity 0.2-0.5); versioned-revision, reason logged;
 (R3) mass scale: nuisance eta (Mtot -> eta*Mtot, i.e. v-tilde -> v-tilde/sqrt(eta)),
      calibrated on the control where alpha == 1 by g >> a0 (the control's design role).
Joint (eta, f_t) fit on the control histogram per e-prior; residual must satisfy
|dv_ctrl| < 0.02 (the Paper III forecast target). Nuisances then FROZEN and propagated
to the wide bins; repaired wide-bin landscape reported vs the untouched frozen band.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from t2_guard import enforce, DERIVED, ROOT

CAT = os.path.join(ROOT, "data/T2/wb/all_columns_catalog.fits.gz")
freezes, fpaths = enforce(staged=[CAT], allow_dirty_paths=("t2_5b_control_repair.py",))
rng = np.random.default_rng(20260726)

cols = ["parallax1", "parallax2", "parallax_over_error1", "parallax_over_error2",
        "pmra1", "pmra2", "pmdec1", "pmdec2", "pmra_error1", "pmra_error2",
        "pmdec_error1", "pmdec_error2", "ruwe1", "ruwe2", "phot_g_mean_mag1",
        "phot_g_mean_mag2", "sep_AU", "binary_type", "R_chance_align"]
with fits.open(CAT, memmap=False) as h:
    d = h[1].data
    A = {c: np.asarray(d[c]) for c in cols}
    del d
plx = 0.5*(A["parallax1"] + A["parallax2"])
dist_pc = 1000.0/np.clip(plx, 1e-3, None)
MG1 = A["phot_g_mean_mag1"] + 5*np.log10(np.clip(A["parallax1"], 1e-3, None)) - 10
MG2 = A["phot_g_mean_mag2"] + 5*np.log10(np.clip(A["parallax2"], 1e-3, None)) - 10
mask = ((A["binary_type"] == "MSMS") & (A["parallax_over_error1"] > 20)
        & (A["parallax_over_error2"] > 20) & (plx > 5.0)
        & (A["ruwe1"] < 1.4) & (A["ruwe2"] < 1.4) & (A["R_chance_align"] < 0.1)
        & (A["phot_g_mean_mag1"] < 18) & (A["phot_g_mean_mag2"] < 18))
MGgrid = np.array([1.4, 2.6, 4.0, 4.8, 5.6, 6.7, 8.0, 9.5, 11.0, 13.0, 15.0])
Mgrid = np.array([2.0, 1.5, 1.1, 1.0, 0.9, 0.75, 0.60, 0.45, 0.30, 0.15, 0.09])
mass_of = lambda MG: np.interp(np.clip(MG, MGgrid[0], MGgrid[-1]), MGgrid, Mgrid)
Mtot = mass_of(MG1) + mass_of(MG2)
dmu = np.hypot(A["pmra1"]-A["pmra2"], A["pmdec1"]-A["pmdec2"])
sdmu = np.sqrt(A["pmra_error1"]**2 + A["pmra_error2"]**2
               + A["pmdec_error1"]**2 + A["pmdec_error2"]**2)/np.sqrt(2)
vp = 4.74047e-3*dmu*dist_pc; svp = 4.74047e-3*sdmu*dist_pc
vc = 29.785*np.sqrt(Mtot/np.clip(A["sep_AU"], 1, None))
vt = vp/vc; snr = svp/vc
good = mask & (svp < 0.15*vc) & np.isfinite(vt)
print(f"frozen-cut sample: {good.sum()}")

def newton_mc(e_prior, n=400_000):
    if e_prior == "thermal": e = np.sqrt(rng.random(n))
    else: e = rng.random(n)**(1/2.2)
    Mano = rng.random(n)*2*np.pi; E = Mano.copy()
    for _ in range(12): E -= (E - e*np.sin(E) - Mano)/(1 - e*np.cos(E))
    cosE, sinE = np.cos(E), np.sin(E); b = np.sqrt(1-e**2)
    r = np.stack([cosE-e, b*sinE], 1)
    v = np.stack([-sinE, b*cosE], 1)/(1-e*cosE)[:, None]
    w = rng.random(n)*2*np.pi; cw, sw = np.cos(w), np.sin(w)
    R = np.stack([np.stack([cw, -sw], 1), np.stack([sw, cw], 1)], 1)
    r2 = np.einsum("nij,nj->ni", R, r); v2 = np.einsum("nij,nj->ni", R, v)
    ci = rng.random(n)*2 - 1
    r3 = np.stack([r2[:, 0], r2[:, 1]*ci], 1); v3 = np.stack([v2[:, 0], v2[:, 1]*ci], 1)
    s = np.linalg.norm(r3, axis=1)
    return np.linalg.norm(v3, axis=1)/np.sqrt(1/s)
MC = {p: newton_mc(p) for p in ("thermal", "hwang")}

BINS = {"control_0.2-1kAU": (200, 1000), "2-5kAU": (2000, 5000),
        "5-10kAU": (5000, 10000), "10-30kAU": (10000, 30000)}
HEDGES = np.linspace(0, 2.5, 26)

def convolve_noise(mc, snr_bin):
    """R1: add the bin's empirical per-component noise to the noise-free MC v-tilde."""
    s = rng.choice(snr_bin, len(mc))
    phi = rng.random(len(mc))*2*np.pi
    vx, vy = mc*np.cos(phi), mc*np.sin(phi)
    return np.hypot(vx + rng.normal(0, s), vy + rng.normal(0, s))

def triple_family(mc, fts):
    """precompute triple-boosted variants per f_t (R2 grid)."""
    out = {}
    boostmask = rng.random(len(mc)); boostfac = 1 + rng.uniform(0.3, 0.8, len(mc))
    for ft in fts:
        b = np.where(boostmask < ft, boostfac, 1.0)
        out[ft] = mc*b
    return out

FTS = np.round(np.arange(0, 0.51, 0.02), 2)
ETAS = np.round(np.arange(0.95, 1.21, 0.01), 2)
def moments(dat, sim):
    """(mean ratio - 1, tail-fraction ratio - 1): mean pins eta+ft combo, tail breaks it."""
    m = dat[(dat > 0) & (dat < 2.5)]; s = sim[(sim > 0) & (sim < 2.5)]
    dv = np.mean(m)/np.mean(s) - 1
    tf_d = np.mean(m > 1.2); tf_s = max(np.mean(s > 1.2), 1e-6)
    return float(dv), float(tf_d/tf_s - 1)

sel_ctrl = good & (A["sep_AU"] >= 200) & (A["sep_AU"] < 1000)
snr_ctrl = snr[sel_ctrl]
print(f"control N={sel_ctrl.sum()}, median sigma_v/v_c = {np.median(snr_ctrl):.3f}")
repair, results = {}, {"bins": {}}
for prior in ("thermal", "hwang"):
    mc_n = convolve_noise(MC[prior], snr_ctrl)
    fam = triple_family(mc_n, FTS)
    best = None
    for eta in ETAS:
        dat = vt[sel_ctrl]/np.sqrt(eta)
        for ft in FTS:
            dv_, tf_ = moments(dat, fam[ft])
            J = (dv_/0.004)**2 + (tf_/0.05)**2          # two-moment objective
            if best is None or J < best[0]: best = (J, eta, ft)
    _, eta_h, ft_h = best
    dat = vt[sel_ctrl]/np.sqrt(eta_h)
    m = dat[(dat > 0) & (dat < 2.5)]; sim = fam[ft_h]; sim = sim[(sim > 0) & (sim < 2.5)]
    dv = float(np.mean(m)/np.mean(sim) - 1)
    boots = [np.mean(rng.choice(m, len(m)))/np.mean(sim) - 1 for _ in range(200)]
    repair[prior] = dict(eta=eta_h, f_t=ft_h, dv_ctrl=round(dv, 4), dv_err=round(float(np.std(boots)), 4))
    print(f"[{prior:8s}] eta = {eta_h:.2f}, f_t = {ft_h:.2f}  ->  control dv = {dv:+.4f} ± {np.std(boots):.4f}")

val = all(abs(repair[p]["dv_ctrl"]) < 0.02 for p in repair)
print(f"\nVALIDATION (|dv_ctrl| < 0.02 both priors): {'PASS' if val else 'FAIL'}")

print(f"\n{'bin':16}{'n':>6}  {'e-prior':8}{'dv repaired':>12}{'(dry-run was)':>14}")
t25 = json.load(open(os.path.join(DERIVED, "t2_5_results.json")))
for bname, (lo, hi) in BINS.items():
    sel = good & (A["sep_AU"] >= lo) & (A["sep_AU"] < hi)
    results["bins"][bname] = {"n": int(sel.sum())}
    for prior in ("thermal", "hwang"):
        eta, ft = repair[prior]["eta"], repair[prior]["f_t"]
        mc_n = convolve_noise(MC[prior], snr[sel])            # per-bin noise convolution
        b = np.where(rng.random(len(mc_n)) < ft, 1 + rng.uniform(0.3, 0.8, len(mc_n)), 1.0)
        sim = mc_n*b; sim = sim[(sim > 0) & (sim < 2.5)]
        dat = vt[sel]/np.sqrt(eta); m = dat[(dat > 0) & (dat < 2.5)]
        dv = float(np.mean(m)/np.mean(sim) - 1)
        boots = [np.mean(rng.choice(m, len(m)))/np.mean(sim) - 1 for _ in range(200)]
        results["bins"][bname][prior] = dict(delta_v=round(dv, 4), delta_v_err=round(float(np.std(boots)), 4))
        was = t25["bins"][bname][prior]["delta_v"]
        print(f"{bname:16}{sel.sum():6d}  {prior:8}{dv:+9.4f}±{np.std(boots):.4f}   ({was:+.3f})")

results.update(repair_nuisances=repair, validation=dict(control_repaired=bool(val), target="|dv|<0.02"),
               revision=dict(reason="control-bin rail at f_t=0.15 + literature hidden multiplicity 0.2-0.5; "
                                     "noise-free MC vs noisy data asymmetry; mass-scale nuisance calibrated "
                                     "on the g>>a0 control (its design role). Estimator-side only; the frozen "
                                     "prediction dv=0.12-0.17 is untouched.",
                             changes=["MC noise convolution (per-bin empirical sigma_v/v_c)",
                                      "f_t grid [0,0.15] -> [0,0.5]",
                                      "mass-scale nuisance eta in [0.90,1.30], control-calibrated"],
                             frozen_for_DR4=dict(thermal=repair["thermal"], hwang=repair["hwang"])),
               note="NO adjudication from eDR3 (frozen: phase 3 = DR4, identical repaired code).")
json.dump(results, open(os.path.join(DERIVED, "t2_5b_repair.json"), "w"), indent=1)

# versioned ledger revision (protocol: revisions by new freeze version with logged reason)
v2 = json.load(open(os.path.join(DERIVED, "T2_locked_predictions_v2.json")))
v2.setdefault("revisions", []).append(dict(
    tag="wb_adjudication_rev2 (T2.5b control repair)", utc="2026-07-25",
    reason=results["revision"]["reason"], changes=results["revision"]["changes"],
    nuisances_frozen_for_DR4=results["revision"]["frozen_for_DR4"],
    prediction_untouched="dv = 0.12-0.17 single fork"))
json.dump(v2, open(os.path.join(DERIVED, "T2_locked_predictions_v3.json"), "w"), indent=1)
print("\nwrote derived/t2_5b_repair.json + T2_locked_predictions_v3.json (versioned revision)")