#!/usr/bin/env python
"""
D19 — the screen-mode boundary problem: deriving C_T = 10.94.

Structural theorems extracted from the matching requirements (each derived, not assumed):
 (i)  omega_gap^2 ∝ X  =>  the confined quantum is a PHONON (gap ∝ c_s);
      particle confinement (omega ∝ hbar/m d^2, X-independent) is EXCLUDED.
 (ii) C_T universal (D18 cancellation) => the gap-setting length is c_s-independent:
      d ∝ n^{-1/3}. Healing- or thermal-length screens break universality (shown below).
 (iii) EXACT SIMPLIFICATION: (2 sqrt2)^{2/3} = 2  =>  C_T = 4 * 3^{1/3} * zeta(3/2)^{2/3};
      in interparticle units zeta cancels: required thickness a = d n^{1/3} with
      a^2 = pi/(4 * 3^{1/3}) = 0.5446 (a = 0.7380) — the screen is athermal pure geometry.
 (iv) Two-scale structure forced: foliation SPACING = l_Pl,ac ∝ xi (entropic side; needed
      for the D18 cancellation via N_A) while profile THICKNESS ∝ n^{-1/3} (matter side).
      Geometric consistency: thickness < spacing in all environments (checked).

Then: the candidate parameter-free profiles vs the exact target
      k_1 = 2 sqrt(pi) 3^{1/6} n^{1/3} = 4.2573 n^{1/3}.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t2_guard import DERIVED, ROOT

HBAR, KB = 1.054571817e-34, 1.380649e-23
EV = 1.78266192e-36
m = 8.5 * EV
ZETA = 2.6124
C_T_REQ = 4 * 3 ** (1 / 3) * ZETA ** (2 / 3)

# --- (iii) verify the exact identity ---
lhs = 6 * (2 * np.sqrt(2) / 3) ** (2 / 3) * ZETA ** (2 / 3)
out = {"identity_check": dict(D18_form=float(lhs), simplified=float(C_T_REQ),
                              statement="6(2sqrt2/3)^{2/3} = 4*3^{1/3} exactly "
                                        "((2sqrt2)^{2/3}=2)")}

# --- (ii) exclusion tests: C_T under alternative thickness rules across environments ---
def Tc_freq(n): return (2 * np.pi / ZETA ** (2 / 3)) * HBAR * n ** (2 / 3) / m
envs = {"galactic": (4.5e13, 1.1e5), "dwarf": (1e12, 1.3e4), "cosmic": (1.45e8, 0.33)}
rules = {}
for rule in ("d=xi", "d=n^-1/3", "d=0.738 n^-1/3"):
    vals = {}
    for k, (n, cs) in envs.items():
        X = m * cs ** 2 / (2 * HBAR)
        if rule == "d=xi":
            d = HBAR / (m * cs)
        elif rule == "d=n^-1/3":
            d = n ** (-1 / 3)
        else:
            d = 0.7380 * n ** (-1 / 3)
        w2 = cs ** 2 * (np.pi / d) ** 2
        vals[k] = float(w2 / (X * Tc_freq(n)))
    rules[rule] = vals
out["exclusion_tests"] = dict(
    results=rules,
    verdict="d=xi: C_T varies by ~8 orders across environments — EXCLUDED (breaks D18 "
            "universality); d ∝ n^{-1/3}: constant — the screen thickness is density-set. "
            "The required a = 0.7380 reproduces C_T = 10.94 in every environment.")

# --- (iv) two-scale geometry check ---
geom = {}
for k, (n, cs) in envs.items():
    xi = HBAR / (m * cs)
    geom[k] = dict(spacing_lPl_ac=float(np.sqrt(12 * np.pi) * xi),
                   thickness=float(0.738 * n ** (-1 / 3)),
                   ratio=float(0.738 * n ** (-1 / 3) / (np.sqrt(12 * np.pi) * xi)))
out["two_scale_geometry"] = dict(values=geom,
                                 note="thickness < spacing everywhere sampled — thin, "
                                      "well-separated screens: geometrically consistent")

# --- candidate parameter-free profiles vs the target ---
kD = (6 * np.pi ** 2) ** (1 / 3)                       # Debye, single branch
cands = {
    "hard-wall d = n^-1/3 (DD)": np.pi * ZETA ** (2 / 3),
    "half-thermal d = lambda_c/2 (DD)": 4 * np.pi,
    "Debye saturation k1 = k_D": kD ** 2 * ZETA ** (2 / 3) / np.pi,
    "k1^3 = 24 pi n (flagged, underived)": (24 * np.pi) ** (2 / 3) * ZETA ** (2 / 3) / np.pi,
}
out["candidates"] = {k: dict(C_T=float(v), dev_pct=float(100 * (v / C_T_REQ - 1)))
                     for k, v in cands.items()}
out["required"] = dict(C_T=float(C_T_REQ), a=0.7380, a2_exact="pi/(4*3^{1/3})",
                       k1_over_n13=float(2 * np.sqrt(np.pi) * 3 ** (1 / 6)),
                       k1_vs_debye=float(2 * np.sqrt(np.pi) * 3 ** (1 / 6) / kD))

out["verdict"] = (
    "NOT CLOSED — REDUCED TO FINAL FORM, with four structural theorems derived en route: "
    "(i) the confined quantum is a phonon (X-linearity); (ii) the screen thickness is "
    "density-set, d ∝ n^{-1/3} (healing/thermal screens excluded by 8-dex universality "
    "violation); (iii) the screen is athermal pure geometry (zeta cancels; "
    "C_T = 4*3^{1/3} zeta^{2/3} exactly); (iv) two-scale structure (entropic spacing l_Pl,ac, "
    "matter thickness ~0.74/n^{1/3}) — geometrically consistent. The boundary problem is now "
    "ONE 1D eigenvalue with its target known to four digits: k_1 = 2 sqrt(pi) 3^{1/6} n^{1/3} "
    "= 4.2573 n^{1/3} (d_eff = 0.7380 interparticle spacings). No canonical parameter-free "
    "profile lands exactly: hard-wall interparticle slab -46%; Debye saturation -16% (the "
    "leading mechanistic candidate: screens where the phonon mode budget saturates); "
    "half-thermal +15%; the 24-pi coincidence (-1.6%) is flagged as underived numerology. "
    "Kill-scope PASSED (a = O(1); an a >> 1 or << 1 would have falsified the screen picture). "
    "What closes it: the self-consistent screen density profile from the entropic sector's "
    "back-reaction (a dark-soliton-like depletion notch) and its lowest Bogoliubov eigenvalue "
    "— a well-posed GPE boundary computation, target 0.7380.")
json.dump(out, open(os.path.join(DERIVED, "b2_screen_eigenvalue.json"), "w"), indent=1)

v2p = os.path.join(DERIVED, "T2_locked_predictions_v2.json")
v2 = json.load(open(v2p))
v2["branches"]["B2"]["screen_eigenvalue"] = dict(
    status="reduced to one 1D eigenvalue: d_eff = 0.7380 n^{-1/3} (k1 = 4.2573 n^{1/3})",
    theorems="phonon confinement; density-set thickness; athermal geometry; two-scale structure",
    C_T_simplified="4*3^{1/3} zeta(3/2)^{2/3}",
    leading_candidate="Debye-saturation boundary (-16%)",
    ref="docs/DM_DERIVATIONS.md D19")
v2["provenance"]["revision_reason"] += " | D19: screen eigenvalue reduced (same cycle)."
json.dump(v2, open(v2p, "w"), indent=1)

print("=== D19: screen-mode boundary problem ===")
print(f"identity: 6(2sqrt2/3)^(2/3) zeta^(2/3) = {lhs:.4f} = 4*3^(1/3) zeta^(2/3) = {C_T_REQ:.4f}")
print("\nexclusion tests (C_T across environments):")
for rule, vals in rules.items():
    print(f"  {rule:18} " + "  ".join(f"{k}={v:.3g}" for k, v in vals.items()))
print("\ncandidate profiles:")
for k, v in out["candidates"].items():
    print(f"  {k:38} C_T = {v['C_T']:7.3f}  ({v['dev_pct']:+.1f}%)")
print(f"\nrequired: C_T = {C_T_REQ:.4f};  a = 0.7380 (a^2 = pi/(4*3^{{1/3}}));  "
      f"k1 = {2*np.sqrt(np.pi)*3**(1/6):.4f} n^(1/3)  ({out['required']['k1_vs_debye']:.3f} x k_D)")
print("\nsaved derived/b2_screen_eigenvalue.json + freeze v2 updated")
