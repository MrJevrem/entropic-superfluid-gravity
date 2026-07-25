#!/usr/bin/env python
"""
D18 — the CW coefficient matching: K_CW = K_locked (the last derivational gap).

Structure:
 1. K_locked (galactic, operational): K = (2 hbar n/3) sqrt(2 hbar/(m c_s^2)), c_s = v_f/sqrt2.
 2. NEGATIVE LEMMA: bulk phonon loops give the LHY X^{5/2} term, not X^{3/2}
    (ratio to P_ent = 1/(5 pi^2 n xi^3) — identical object to D13.1's eps_loop):
    the stiffness source is IRREDUCIBLE to the condensate (third irreducibility result).
 3. Screen-CW matching: K_CW = hbar gamma_c^{3/2} N_A/(12 pi c_s^2), c_bar = c_s (forced),
    N_A = 1/l_Pl,ac = 1/(sqrt(12 pi) xi) (D9's own cutoff). Invert for gamma_c.
 4. THE DISCOVERY: the c_s-dependence cancels IDENTICALLY in the inversion:
       gamma_c = A^{2/3} * hbar n^{2/3}/m,  A = (12 pi)^{3/2} * 2 sqrt2/3
    i.e. gamma_c is proportional to the condensate's OWN critical-temperature frequency:
       hbar gamma_c = C_T * k_B T_c(n),  C_T = 6 (2 sqrt2/3)^{2/3} zeta(3/2)^{2/3} = 10.94
    — an environment-independent pure number built from D9's 12 pi and BEC statistics.
    K_CW = K_locked holds in EVERY environment simultaneously iff the screen-mode gap law is
       omega_gap^2 = C_T (k_B T_c/hbar) X.
 5. Kill-scope: gamma_c could have demanded an unavailable scale (>> Compton or << chemical).
    It lands at ~11 T_c — the condensate's coherence scale: physically natural.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t2_guard import DERIVED, ROOT

HBAR, KB = 1.054571817e-34, 1.380649e-23
EV = 1.78266192e-36
m = 8.5 * EV
ZETA = 2.612

def K_locked(n, cs):
    return (2 * HBAR * n / 3) * np.sqrt(2 * HBAR / (m * cs ** 2))

def gamma_c_required(n, cs):
    K = K_locked(n, cs)
    xi = HBAR / (m * cs)
    N_A = 1 / (np.sqrt(12 * np.pi) * xi)
    return (K * 12 * np.pi * cs ** 2 / (HBAR * N_A)) ** (2 / 3)

def Tc_freq(n):                                   # k_B T_c / hbar
    return (2 * np.pi / ZETA ** (2 / 3)) * HBAR * n ** (2 / 3) / m

# fiducial + environment sweep (the cancellation test)
envs = {"galactic_mid": (4.5e13, 1.1e5), "outskirts": (6.6e11, 1.1e5 * 0.6),
        "dwarf": (1e12, 1.3e4), "cosmic": (1.45e8, 0.33)}
sweep = {}
for k, (n, cs) in envs.items():
    g = gamma_c_required(n, cs)
    sweep[k] = dict(gamma_c=float(g), ratio_to_Tc=float(g / Tc_freq(n)))
C_T_exact = 6 * (2 * np.sqrt(2) / 3) ** (2 / 3) * ZETA ** (2 / 3)
K_fid = K_locked(4.5e13, 1.1e5)
g_fid = gamma_c_required(4.5e13, 1.1e5)

# negative lemma cross-check: LHY/P_ent vs D13.1 eps_loop
n_g, cs_g = 4.5e13, 1.1e5
xi_g = HBAR / (m * cs_g)
lhy_ratio = 1 / (5 * np.pi ** 2 * n_g * xi_g ** 3)
out = dict(
    K_locked_fiducial=float(K_fid),
    K_band="x[0.1,10] from (n_gal, v_flat) choices",
    negative_lemma=dict(
        statement="bulk phonon loops -> LHY term ∝ c_s^5 ∝ X^{5/2}, NOT X^{3/2}: the "
                  "stiffness source is irreducible to the condensate (third irreducibility "
                  "result, joining D11/D12's noise irreducibility)",
        lhy_over_P_ent=float(lhy_ratio), d13_eps_loop=3.55e-4,
        cross_check="same object within O(1) conventions — D13.1's loop parameter IS the "
                    "bulk X^{5/2} correction"),
    inversion=dict(gamma_c_fiducial=float(g_fid), hbar_gamma_c_meV=float(HBAR * g_fid / EV / 1e-3 * 8.5 / 8.5 / 1e0) if False else float(HBAR * g_fid / 1.602e-19 * 1e3),
                   compton_freq=float(m * 9e16 / HBAR), chemical_freq=float(m * cs_g ** 2 / HBAR / 2) * 2),
    cancellation=dict(
        sweep=sweep,
        statement="gamma_c/(k_B T_c/hbar) is the SAME constant in every environment — the "
                  "c_s-dependence cancels identically in the inversion",
        C_T_exact_formula="6 (2 sqrt2/3)^{2/3} zeta(3/2)^{2/3}",
        C_T=float(C_T_exact)),
    gap_law="omega_gap^2 = C_T (k_B T_c/hbar) X   [frozen]",
    verdict=(
        "MATCHING CLOSES CONDITIONALLY — AND THE CONDITION IS UNIVERSAL. The inversion's "
        "c_s-dependence cancels exactly, so K_CW = K_locked holds in every environment "
        "simultaneously iff the screen-mode gap rides the condensate's critical temperature "
        f"with the universal coefficient C_T = {C_T_exact:.2f} — a pure number composed of "
        "D9's counted 12pi and BEC statistics (zeta(3/2)). hbar*gamma_c ~ 11 k_B T_c: "
        "physically natural (the screens' gap = the coherence scale), mid-spectrum "
        "(4e4 below Compton, 180x above chemical) — the kill-scope test passes. "
        "Residual [O] shrinks to ONE pure number: derive C_T = 10.94 from the screen-mode "
        "boundary problem. Emergence status: form (D17, CW-natural), environment-consistency "
        "(D18, exact cancellation), coefficient (one dimensionless number away)."))
json.dump(out, open(os.path.join(DERIVED, "b2_cw_matching.json"), "w"), indent=1)

v2p = os.path.join(DERIVED, "T2_locked_predictions_v2.json")
v2 = json.load(open(v2p))
v2["branches"]["B2"]["cw_matching"] = dict(
    status="closes conditionally: gap law omega_gap^2 = C_T (k_B T_c/hbar) X frozen",
    C_T=round(C_T_exact, 2), C_T_formula="6(2sqrt2/3)^{2/3} zeta(3/2)^{2/3}",
    negative_lemma="bulk gives X^{5/2} (LHY) — third irreducibility",
    residual="derive C_T from screen boundary problem",
    ref="docs/DM_DERIVATIONS.md D18")
v2["provenance"]["revision_reason"] += " | D18: CW matching -> universal gap law (same cycle)."
json.dump(v2, open(v2p, "w"), indent=1)

print("=== D18: CW coefficient matching ===")
print(f"K_locked (fiducial) = {K_fid:.3e} SI")
print(f"negative lemma: LHY/P_ent = {lhy_ratio:.2e} vs D13.1 eps = 3.6e-4 (same object, O(1))")
print(f"required gamma_c = {g_fid:.3e} s^-1  (hbar*gamma_c = {HBAR*g_fid/1.602e-19*1e3:.2f} meV)")
print(f"  vs Compton {m*9e16/HBAR:.1e} | chemical {m*cs_g**2/HBAR:.1e}  -> mid-spectrum OK")
print("cancellation test (gamma_c / (k_B T_c/hbar) across environments):")
for k, v in sweep.items():
    print(f"  {k:14} gamma_c = {v['gamma_c']:.3e}   ratio to T_c-freq = {v['ratio_to_Tc']:.4f}")
print(f"C_T exact = 6(2√2/3)^(2/3) ζ(3/2)^(2/3) = {C_T_exact:.4f}")
print("\n=> gap law omega_gap^2 = C_T (k_B T_c/hbar) X  closes K_CW = K_locked universally")
print("saved derived/b2_cw_matching.json + freeze v2 updated")
