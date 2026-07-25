#!/usr/bin/env python
"""
T3.0 — THE FREEZE for the sigma-resolved RAR test (T3_PLAN).

Everything below is fixed BEFORE any Tier-1..7 dataset is downloaded. The
sigma_crit kernel is NOT re-derived: it is reconstructed from the locked D5
entry (b2_derived.json / T2_locked_predictions_v2.json) and verified to
reproduce the locked fiducial 600 km/s at m = 8.453 eV, rho_ref = 1e-22 kg/m^3.

REVISION NOTE (logged here, plan not silently edited): T3_PLAN quoted a search
window [400, 900] km/s and an elliptical-shelf kill "below 400" — those numbers
came from the papers' narrative "500-800" band, which the ledger shows is the
DENSITY-CONVENTION spread at FIXED fiducial m (rho in [0.5, 2] x rho_ref =>
476-756 km/s), not the m-window mapping. The locked kernel maps the full
m-window [5.76, 16.30] eV to sigma_crit in [250, 1000] km/s at rho_ref.
This freeze carries the kernel-correct windows; the plan discrepancy is
documented in ledger["revision_note"].
"""
import json, os, subprocess, sys
from datetime import datetime, timezone
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import DERIVED, ROOT

H_PLANCK = 6.62607015e-34
ZETA32 = 2.6123753486854883
EV = 1.602176634e-19
C = 2.99792458e8

# ---- reconstruct the locked D5 kernel and verify ----
d5 = json.load(open(os.path.join(DERIVED, "b2_derived.json")))["derivations"]["D5_phase_boundary"]
RHO_REF = d5["rho_halo"]                       # 1e-22 kg/m^3 (locked)
M_FID, M_MIN, M_MAX = (d5["m_eV"][k] for k in ("fiducial", "min", "max"))
SIG_FID_LOCKED = d5["sigma_crit_fiducial_km_s"]  # 600.0 (locked)

def sigma_crit_kms(m_eV, rho=RHO_REF):
    """Locked kernel: n lambda^3 = zeta(3/2), lambda = h/(m sigma), n = rho/m."""
    m = m_eV * EV / C**2
    return (rho * H_PLANCK**3 / (ZETA32 * m**4)) ** (1.0 / 3.0) / 1e3

rec = sigma_crit_kms(M_FID)
assert abs(rec - SIG_FID_LOCKED) < 2.0, f"kernel reconstruction failed: {rec} vs {SIG_FID_LOCKED}"

RHO_BAND = (0.5, 2.0)   # frozen density-convention fork (explains the papers' 500-800 band)
sig_lo_m, sig_hi_m = sigma_crit_kms(M_MAX), sigma_crit_kms(M_MIN)      # m-window at rho_ref
band_fid = tuple(round(sigma_crit_kms(M_FID, f * RHO_REF), 1) for f in RHO_BAND)
outer = (round(sig_lo_m * RHO_BAND[0] ** (1 / 3), 0), round(sig_hi_m * RHO_BAND[1] ** (1 / 3), 0))

A0_LOCKED = 1.0801e-10   # m/s^2, c^2 sqrt(Lambda/24pi), Planck 2018 (N2)

git = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True)
ledger = {
  "provenance": {"git_commit": git.stdout.strip(),
                 "utc": datetime.now(timezone.utc).isoformat(),
                 "script": "analysis/t3_predictions_freeze.py",
                 "plan": "docs/T3_PLAN_sigma_resolved_RAR.md"},
  "revision_note": ("Plan windows [400,900] and shelf-kill '<400' superseded at freeze: the "
                    "narrative 500-800 band is the rho-convention spread at fixed m_fid "
                    "(476-756 km/s), NOT the m-window mapping. Kernel-correct values below."),
  "kernel": {"form": "sigma_crit = [rho h^3 / (zeta(3/2) m^4)]^(1/3), lambda_dB = h/(m sigma), n = rho/m",
             "rho_ref_kg_m3": RHO_REF, "rho_convention_fork": RHO_BAND,
             "m_eV": {"min": M_MIN, "fiducial": M_FID, "max": M_MAX},
             "sigma_crit_km_s": {"at_m_fid": round(rec, 1),
                                 "m_window_at_rho_ref": [round(sig_lo_m, 0), round(sig_hi_m, 0)],
                                 "rho_band_at_m_fid": band_fid,
                                 "outer_envelope": outer},
             "inversion": "m/m_fid = (sigma_b/600 km/s)^(-3/4) * (rho/rho_ref)^(1/4)",
             "inversion_note": "rho nuisance enters at 1/4 power: factor-2 rho error -> 19% in m"},
  "statistic": {
      "D_def": "mean over radial points of log10(g_obs) - log10(g_RAR(g_bar)), per system",
      "g_bar_window_m_s2": [1e-12, 1e-10], "min_points_in_window": 3,
      "g_RAR_primary": {"form": "g_bar / (1 - exp(-sqrt(g_bar/a0)))", "a0_m_s2": A0_LOCKED,
                        "tag": "locked theory value [N2]"},
      "g_RAR_fork_empirical": {"a0_m_s2": 1.20e-10, "tag": "canonical g-dagger (comparator fork)"},
      "cosmology": {"H0": 67.36, "Om": 0.315, "flat": True}},
  "sigma_proxies": {
      "members": "1D line-of-sight dispersion, gapper estimator, >=10 members",
      "TX": {"form": "sigma = 400.4 * sqrt(T_X/keV) km/s  (beta_spec = 1)",
             "fork": "beta_spec = 0.9 -> multiply by 0.949"},
      "lensing_mass": {"form": "sigma = 1082.9 * (E(z) M200c / 1e15 Msun)^0.336 km/s (Evrard+08)"},
      "rule": "method-matched comparisons primary; cross-method overlap reported as systematic"},
  "predictions": {
      "P1_condensed_floor": {"claim": "|mean D| < 0.10 dex for relaxed systems with sigma < 200 km/s "
                                       "(guaranteed-condensed for ALL m in window x rho band)",
                             "note": "the 250-400 km/s shelf is MEASUREMENT territory, not a kill zone: "
                                      "a break there implies high-m carrier via the inversion"},
      "P2_normal_excess": {"claim": "mean D > 0.15 dex for sigma > 1000 km/s; strict sub-check sigma > 1260 "
                                     "(guaranteed-normal beyond the m_min x rho-high corner)"},
      "P3_break": {"claim": "single transition in D(sigma) with center sigma_b in the outer envelope "
                             f"{list(outer)} km/s; fiducial forecast {round(rec,0)} km/s "
                             f"(rho band {list(band_fid)})",
                   "deliverable": "posterior on sigma_b -> posterior on m via the inversion"},
      "P4_discriminator": {"claim": "dD/dln(sigma) > 0 at >=3 sigma in straddle bins AT FIXED M; "
                                     "dD/dln(M) at fixed sigma consistent with 0 (<2 sigma); "
                                     "sigma_b independent of z within errors"}},
  "kill_conditions": [
      "P4 inverted: D tracks M at fixed sigma at >=3 sigma significance",
      f"break center sigma_b outside {list(outer)} km/s at >=2 sigma",
      "no transition detected with mean D(>1260) - D(<200) < 0.10 dex (P2 fails)",
      "systematic RAR failure |mean D| > 0.15 dex in relaxed systems below 200 km/s (P1 fails)"],
  "forks_frozen": {
      "rho_convention": RHO_BAND, "beta_spec": [1.0, 0.9],
      "a0": ["locked 1.0801e-10", "empirical 1.20e-10"],
      "sample": ["relaxed-only (published flags; cool-core proxy where absent)", "all-systems"],
      "reporting_rule": "every result reported under all forks; no post-hoc fork selection"},
  "staging_manifest_order": [
      "X-COP tables", "ACCEPT profile library", "Sun2009+Lovisari group profiles",
      "HeCS caustic profiles", "GAMA G3C", "ETG dynamics tables (shelf)",
      "WL mass calibrations", "(stretch) eRASS1", "(stretch) CLASH/HFF"],
  "guard": "analysis/t3_guard.py — all comparison data must be staged AFTER this freeze"}

out = os.path.join(DERIVED, "T3_locked_predictions.json")
json.dump(ledger, open(out, "w"), indent=1)
print("=== T3.0 FREEZE ===")
print(f"kernel verified: sigma_crit(m_fid) = {rec:.1f} km/s (locked: {SIG_FID_LOCKED})")
print(f"m-window -> sigma_crit [{sig_lo_m:.0f}, {sig_hi_m:.0f}] km/s at rho_ref; "
      f"rho band at m_fid: {band_fid}; outer envelope: {outer}")
print(f"commit {ledger['provenance']['git_commit'][:8]}  utc {ledger['provenance']['utc']}")
print(f"wrote {out}")