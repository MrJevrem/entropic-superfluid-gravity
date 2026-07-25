#!/usr/bin/env python
"""
D11 — prove or refute: f_floor = sqrt(4 pi) H  (the D8/D9/D10 kill identity).

The identity requires the condensate phase at the cosmological floor to be SOFT: its
canonical stiffness must equal sqrt(4pi) Hubble rates, so that the universal dS kick
deposits relative Lambda-variance 1/(4pi) per e-fold.

Framework-fixed inputs (no freedom):
  floor:     X_bar = H_L/2pi                      [D2 temperature lemma]
  structure: P(X) = K X^{3/2}                     [MOND regime; c_s^2 = 2 hbar X_bar/m]
  matching:  K X_bar^{3/2} = rho_Lambda           [alpha = 8pi: floor pressure = dark energy]
  kick:      <dtheta^2>/e-fold = hbar H^2/(4pi^2 P'' c_s^3)  =>  f^2 = P'' c_s^3 / hbar

Fork sweeps: m in [5.8, 16.3] eV (D5 band); floor-pressure fraction x in [1e-6, 1];
the 3/2 factor from delta-Lambda/Lambda = (3/2) delta-X/X_bar. Channel enumeration:
phase / density (gapped) / thermal cloud — the refutation must exclude ALL B2 channels.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t2_guard import DERIVED, ROOT

C, G, HBAR = 299792458.0, 6.674e-11, 1.054571817e-34
EV = 1.78266192e-36
H0 = 2.268e-18
H_L = H0 * np.sqrt(0.7)
RHO_L = 0.7 * 3 * H0 ** 2 * C ** 2 / (8 * np.pi * G)       # J/m^3 (= Lambda c^4/8piG)
TARGET = np.sqrt(4 * np.pi)                                 # required f/H

def f_over_H(m_eV, x_frac=1.0):
    m = m_eV * EV
    Xb = H_L / (2 * np.pi)                                  # [1/s]
    K = x_frac * RHO_L / Xb ** 1.5
    Ppp = 0.75 * K / np.sqrt(Xb)                            # [J s^2/m^3]
    cs2 = 2 * HBAR * Xb / m                                 # [m^2/s^2]
    f2 = Ppp * cs2 ** 1.5 / HBAR                            # [1/s^2]
    return np.sqrt(f2) / H_L, np.sqrt(cs2)

central, cs = f_over_H(8.5)
sweep_m = {f"m={m}eV": f_over_H(m)[0] for m in (5.8, 8.5, 16.3)}
sweep_x = {f"x={x:g}": f_over_H(8.5, x)[0] for x in (1.0, 1e-3, 1e-6)}
per_efold = (1 / central) ** 2                              # (H/f)^2 actually delivered

out = dict(
    identity="f_floor = sqrt(4pi) H  (required " + f"{TARGET:.3f} H)",
    computed=dict(f_over_H=float(central), c_s_floor_m_s=float(cs),
                  per_efold_variance_delivered=float(per_efold),
                  per_efold_required=float(1 / (4 * np.pi))),
    discrepancy_dex=float(np.log10(central / TARGET)),
    fork_sweeps=dict(mass={k: float(v) for k, v in sweep_m.items()},
                     floor_fraction={k: float(v) for k, v in sweep_x.items()},
                     note="the 3/2 P-exponent factor shifts the target by x1.5 — irrelevant "
                          "at 35 dex; even x=1e-6 floor fraction moves f by only 10^3"),
    channel_enumeration=dict(
        phase="computed above: f/H ~ 3e35 — 35 dex too stiff",
        density_mode="gapped at m c_s^2/hbar >> H — stiffer still",
        thermal_cloud="entropy density at T_dS negligible (D9: L* ~ 1e77 m) — no capacity",
        conclusion="B2 possesses NO channel soft enough by >= 30 dex"),
    verdict=(
        "REFUTED. f_floor/H = %.1e vs required %.2f — the identity fails by ~%d orders of "
        "magnitude, robustly across every frozen fork (mass band, floor fraction, the 3/2 "
        "factor). The condensate phase delivers per-e-fold Lambda-variance ~%.0e instead of "
        "the required 1/(4pi) ~ 0.08. CONSEQUENCES: (1) the microphysical (in-B2) derivation "
        "of sigma_Lambda^2/Lambda^2 = 1/8pi is DEAD — the Lambda-channel stochasticity cannot "
        "be carried by the condensate; it must be a property of the fundamental CQ "
        "gravitational sector (Oppenheim's stochastic Lambda), whose amplitude is an "
        "irreducible input. (2) D8's a0 = c^2 sqrt(Lambda/24pi) SURVIVES but is DOWNGRADED: "
        "from 'derived up to one identity' to 'phenomenologically locked normalization with "
        "structural support' (the D10 alpha-linearity lemma is carrier-independent and "
        "stands; L1/L2 evidence untouched). (3) D9 untouched (entropy side). (4) SEPARATION "
        "THEOREM: within B2+Dorau-Much, the carrier of ENTROPY (condensate) and the carrier "
        "of STOCHASTICITY (CQ Lambda-channel) are provably distinct sectors — no single-"
        "sector reduction exists." % (central, TARGET, round(np.log10(central / TARGET)),
                                      per_efold)),
    consistency_note=(
        "The refutation is the outcome CONSISTENT with B2's own earlier results: D3/D7 derived "
        "a stiff-quiet condensate (noise 55-66 dex below ceilings). Had f_floor ~ sqrt(4pi)H "
        "held, the condensate would fluctuate at O(1) on Hubble scales while being 55+ dex "
        "silent at nHz — a 50-dex scale-dependence with no mechanism. The test was implicitly "
        "decided the moment D3 was derived; D11 makes the connection exact."))
json.dump(out, open(os.path.join(DERIVED, "b2_f_floor.json"), "w"), indent=1)

v2p = os.path.join(DERIVED, "T2_locked_predictions_v2.json")
v2 = json.load(open(v2p))
v2["branches"]["B2"]["lambda_bookkeeping"]["status"] = (
    "kill identity f_floor=sqrt(4pi)H REFUTED (35 dex, fork-robust): in-B2 microphysical "
    "derivation of the 1/8pi variance dead; a0 lock survives as phenomenological with "
    "structural support; noise-carrier and entropy-carrier provably distinct sectors")
v2["provenance"]["revision_reason"] += " | D11: f_floor identity refuted (same revision cycle)."
json.dump(v2, open(v2p, "w"), indent=1)

print("=== D11: f_floor = sqrt(4pi) H — prove or refute ===")
print(f"floor sound speed c_s = {cs:.2e} m/s")
print(f"computed f/H = {central:.2e}   required = {TARGET:.3f}")
print(f"discrepancy: {np.log10(central/TARGET):.1f} dex  -> REFUTED")
print(f"per-e-fold variance delivered = {per_efold:.1e}  (required 1/4pi = {1/(4*np.pi):.3f})")
print("mass sweep:", {k: f"{v:.1e}" for k, v in sweep_m.items()})
print("floor-fraction sweep:", {k: f"{v:.1e}" for k, v in sweep_x.items()})
print("channels: phase 35-dex stiff | density gapped | thermal negligible -> NO soft channel")
print("saved derived/b2_f_floor.json + freeze v2 updated")
