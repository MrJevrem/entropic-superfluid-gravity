#!/usr/bin/env python
"""
D21 — the many-body mode-selection computation (the last number).

Target (exact, from the lock): k_1 = 2 sqrt(pi) 3^{1/6} n^{1/3}, i.e. mode budget
nu* = k_1^3/(6 pi^2 n) = 4/sqrt(3 pi) = 1.30293 modes/particle;  C_T = zeta^{2/3} k_1^2/(pi n^{2/3}).

Movement 1 — WEIGHTED-BUDGET EXCLUSION (computed): any spectral weighting (structure-factor
  S(k), f-sum) drags in xi*n^{1/3}, which varies ~4 dex across environments => k_1/n^{1/3}
  drifts => excluded by D18 universality. The selection functional must be pure k-space
  geometry. (The many-body derivation of D19's athermal theorem.)
Movement 2 — COUNTING NO-GO (exact): nu* = 4/sqrt(3 pi) is irrational. THEOREM: no sharp
  mode-count (rational nu: 1, 4/3, 3/2, 2, ...) closes C_T exactly. Sharp Debye (nu = 1)
  remains the leading approximation at -16%. The closure requires a SMEARED branch edge.
Movement 3 — EDGE-PROFILE DETERMINATION: for parameter-free reference edges compute the
  selection; the exact closure pins ONE moment of the edge profile:
  <k^2>^{1/2}/k_g = 1.2011 (Gaussian edge would give 1.2247: +2.0% in k, +3.9% in C_T —
  the best parameter-free estimate). The last number is a MATERIAL CONSTANT of the quintic
  superfluid's branch termination — measurable in a cold-atom analog (mu ∝ n^2 EOS).
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from t2_guard import DERIVED, ROOT

HBAR = 1.054571817e-34
EV = 1.78266192e-36
m = 8.5 * EV
ZETA = 2.6124
K1_REQ = 2 * np.sqrt(np.pi) * 3 ** (1 / 6)          # 4.2572, units n^{1/3}
NU_REQ = K1_REQ ** 3 / (6 * np.pi ** 2)             # 4/sqrt(3 pi)
KD = (6 * np.pi ** 2) ** (1 / 3)
def C_T(k1_over_n13): return ZETA ** (2 / 3) * k1_over_n13 ** 2 / np.pi
C_T_REQ = C_T(K1_REQ)

# ---------------- Movement 1: weighted-budget exclusion ----------------
def k1_S_weighted(xi_n13):
    """solve (1/2pi^2) int_0^{k1} k^2 S_Bog(k) dk = n, S = (k xi/2)/sqrt(1+(k xi/2)^2);
       returns k1/n^{1/3}."""
    target = 2 * np.pi ** 2 * xi_n13 ** 3            # in kappa = k*xi units
    F = lambda kap: quad(lambda x: x ** 3 / (2 * np.sqrt(1 + x ** 2 / 4)), 0, kap)[0] - target
    lo, hi = 1e-3, 10.0
    while F(hi) < 0: hi *= 2
    kap1 = brentq(F, lo, hi)
    return kap1 / xi_n13

envs = {"galactic": 2.25, "dwarf": 19.0, "cosmic": 1.108e4}   # xi * n^{1/3}
m1 = {k: float(k1_S_weighted(v)) for k, v in envs.items()}

# ---------------- Movement 2: counting no-go ----------------
rationals = {"nu=1 (Debye)": 1.0, "nu=4/3": 4 / 3, "nu=3/2": 1.5, "nu=2 (two chiralities)": 2.0}
m2 = {k: dict(k1=float((6 * np.pi ** 2 * v) ** (1 / 3)),
              CT_dev_pct=float(100 * (C_T((6 * np.pi ** 2 * v) ** (1 / 3)) / C_T_REQ - 1)))
      for k, v in rationals.items()}

# ---------------- Movement 3: edge profiles ----------------
kg = (8 * np.pi ** 1.5) ** (1 / 3)                   # Gaussian budget: int e^{-k^2/kg^2} = n
edges = {
    "sharp Debye edge (k1 = k_D)": KD,
    "Gaussian edge, rms selection": np.sqrt(1.5) * kg,
    "Gaussian edge, mean selection": 2 * kg / np.sqrt(np.pi),
    "sharp edge, rms selection": np.sqrt(0.6) * KD,
}
m3 = {k: dict(k1=float(v), k_dev_pct=float(100 * (v / K1_REQ - 1)),
              CT_dev_pct=float(100 * (C_T(v) / C_T_REQ - 1))) for k, v in edges.items()}
moment_req = K1_REQ / kg                              # required <k^2>^{1/2}/k_g

out = dict(
    target=dict(k1_over_n13=float(K1_REQ), nu_star=float(NU_REQ),
                nu_star_exact="4/sqrt(3 pi)", C_T=float(C_T_REQ)),
    movement1=dict(k1_over_n13_by_env=m1,
                   verdict="S(k)-weighted budget drifts from "
                           f"{m1['galactic']:.2f} to {m1['cosmic']:.2f} across environments "
                           "— EXCLUDED by D18 universality: the selection functional must be "
                           "pure k-space geometry (weight = fixed dimensionless profile). "
                           "D19's athermal theorem, derived at the functional level."),
    movement2=dict(candidates=m2,
                   theorem="nu* = 4/sqrt(3 pi) is irrational: NO sharp mode-count closes C_T "
                           "exactly. Sharp Debye stays the leading approximation (-16%); the "
                           "closure requires a smeared branch-termination edge."),
    movement3=dict(edges=m3, gaussian_budget_scale=float(kg),
                   required_moment=float(moment_req), gaussian_moment=float(np.sqrt(1.5)),
                   verdict=f"best parameter-free prescription: Gaussian edge + rms selection "
                           f"= +2.0% in k (+3.9% in C_T). Exact closure pins ONE moment of "
                           f"the edge profile: <k^2>^{{1/2}}/k_g = {moment_req:.4f} "
                           f"(Gaussian: 1.2247)."),
    final_verdict=(
        "THE LAST NUMBER IS A MATERIAL CONSTANT. Three-movement result: (1) weighted budgets "
        "are excluded by universality — the selection is pure k-geometry; (2) no sharp count "
        "closes it — nu* = 4/sqrt(3pi) is irrational, so the branch edge is smeared; (3) the "
        "exact closure pins one moment of the edge profile (1.2011 vs Gaussian 1.2247 — the "
        "parameter-free estimate lands within 2%). The coefficient C_T is thus exactly known "
        "from the lock, mechanically localized to the phonon-branch termination edge of the "
        "quintic superfluid, approximated parameter-free to 2% in k, and NOT derivable by "
        "counting at any order — it is the analog of where helium's phonon-roton branch ends: "
        "a measured property of the quantum fluid. FORWARD PATH (experimental): a cold-atom "
        "condensate with engineered mu ∝ n^2 EOS (three-body-dominated interactions) realizes "
        "the same universality class; its branch-edge moment IS C_T up to statistics factors "
        "— the dark-matter force-law coefficient is measurable in a laboratory analog. The "
        "derivation program D0-D21 ends by converting its final unknown into an experimental "
        "proposal."))
json.dump(out, open(os.path.join(DERIVED, "b2_mode_selection.json"), "w"), indent=1)

v2p = os.path.join(DERIVED, "T2_locked_predictions_v2.json")
v2 = json.load(open(v2p))
v2["branches"]["B2"]["mode_selection"] = dict(
    nu_star="4/sqrt(3pi) = 1.30293 (irrational: counting no-go)",
    weighted_budgets="excluded by universality",
    best_parameter_free="Gaussian edge + rms: +2.0% in k",
    residual="one edge-profile moment = 1.2011; measurable in mu∝n^2 cold-atom analog",
    ref="docs/DM_DERIVATIONS.md D21")
v2["provenance"]["revision_reason"] += " | D21: mode selection — material-constant verdict (same cycle)."
json.dump(v2, open(v2p, "w"), indent=1)

print("=== D21: many-body mode selection ===")
print(f"target: k1 = {K1_REQ:.4f} n^(1/3);  nu* = {NU_REQ:.5f} = 4/sqrt(3pi);  C_T = {C_T_REQ:.4f}")
print("\nMovement 1 — weighted-budget exclusion (k1/n^(1/3) by environment):")
for k, v in m1.items(): print(f"  {k:9} {v:.3f}")
print("\nMovement 2 — sharp-count candidates:")
for k, v in m2.items(): print(f"  {k:24} k1 = {v['k1']:.4f}   C_T dev = {v['CT_dev_pct']:+.1f}%")
print("\nMovement 3 — edge profiles:")
for k, v in m3.items(): print(f"  {k:32} k1 = {v['k1']:.4f}  ({v['k_dev_pct']:+.2f}% k, "
                              f"{v['CT_dev_pct']:+.1f}% C_T)")
print(f"\nrequired edge moment <k^2>^(1/2)/k_g = {moment_req:.4f}  (Gaussian: {np.sqrt(1.5):.4f})")
print("\nsaved derived/b2_mode_selection.json + freeze v2 updated")
