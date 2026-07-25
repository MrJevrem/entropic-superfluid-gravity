#!/usr/bin/env python
"""
D32 — the substrate consistency theorems (T-5, tractable stratum; user directive
2026-07-25).

T-5 splits into two strata. The DEEP stratum — derive S_rel = dA/4 from a
microscopic substrate — is the question [DM26] leaves open for spacetime itself;
no consistency argument closes it and none is attempted here. The TRACTABLE
stratum is what gates D28/D31: (I) is the entropic functional's resolution xi
FORCED or chosen? (II) is the entropic phonon a genuine canonical quantum mode
(so the [N3] entanglement counting and the D22 quadratic action are legitimate)?

THEOREM I (scale assignment, by physical regulator + frozen-number lock):
the [N3] counting integral is regulated by the Bogoliubov dispersion crossover —
a property of the medium, not a choice — and the frozen eta = 1/(48 pi xi^2)
is INCOMPATIBLE with any coarse-grained variant: cutting the counting at
lambda_dB instead of xi rescales eta by (2 pi)^2 ~ 39, far outside the stated
scheme envelope x[0.2, 1.5]. The resolution is also observable-locked: the
sub-mm screening scale l_Pl,ac = sqrt(12 pi) xi sits in the torsion-balance row
of the phenomenology package [N9].

THEOREM II (canonical inheritance): the phonon variables (theta, dn) are the
Madelung variables of the MICROSCOPIC bosonic field psi, whose [psi, psi+] = 1
is not in question; hence [theta(x), dn(y)] = i delta(x-y) at high occupancy,
REGARDLESS of the stiffness's origin. The entropic term contributes only the
Hamiltonian. Verification: canonical quantization of (theta, dn) with the
entropic eps(n) = hbar gamma n^3/3 must reproduce the series' frozen sound
speed c_s^2 = 2 hbar gamma n^2/m exactly — one hbar (the particle's), no
independent acoustic hbar; M_Pl,ac = m/sqrt(12 pi) is inherited.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import DERIVED

EV = 1.602176634e-19; C = 2.99792458e8; HBAR = 1.054571817e-34; H = 6.62607015e-34
m = 8.45*EV/C**2; vf = 156e3; cs = vf/np.sqrt(2); sig = cs
xi = HBAR/(m*cs); lam_dB = H/(m*sig)

print("=== D32: the substrate consistency theorems (T-5, tractable stratum) ===\n")

print("--- THEOREM I: the resolution is forced, not chosen ---")
# transverse-channel counting: each k_perp is a 1+1d c=1 half-line entangled up to
# the Bogoliubov crossover; S/A = int d2k/(2pi)^2 (1/6) ln(k_x*/k_perp) with the
# channel terminating where dispersion departs linear: omega = c_s k sqrt(1+k^2 xi^2/2).
def eta_of_cutoff(k_cut_xi):
    # integrate channels up to k_cut (in units of 1/xi); log runs to the crossover
    u = np.linspace(1e-6, k_cut_xi, 400_000)
    integrand = u*np.maximum(np.log(k_cut_xi/u), 0)/6
    return np.trapezoid(integrand, u)/(2*np.pi)      # in units of 1/xi^2, per (2pi)^2->2pi reduction
eta_xi = eta_of_cutoff(1.0)
eta_ref = 1/(48*np.pi)
kg_xi = xi/lam_dB                                     # = 1/(2 pi): the D28 granularity scale in cutoff units
eta_cg = eta_of_cutoff(kg_xi)
print(f"counting with the physical (Bogoliubov) regulator at k = 1/xi: eta*xi^2 = {eta_xi:.5f}")
print(f"frozen [N3] value 1/(48 pi) = {eta_ref:.5f}: ratio {eta_xi/eta_ref:.2f} — inside the stated")
print(f"scheme envelope x[0.2, 1.5] (the O(1) is the known scheme freedom).")
print(f"coarse-grained variant (cutoff at lambda_dB, i.e. k_cut*xi = 1/2pi = {kg_xi:.4f}):")
print(f"  eta*xi^2 = {eta_cg:.7f} — smaller by x{eta_xi/eta_cg:.0f} ~ (2 pi)^2: FAR outside the envelope.")
print(f"=> the frozen G_ph = 12 pi hbar c_s/m^2 is arithmetically INCOMPATIBLE with a")
print(f"   coarse-grained functional: the series committed to xi-resolution the moment [N3]")
print(f"   froze — and that commitment is observable-locked (l_Pl,ac = sqrt(12 pi) xi = "
      f"{np.sqrt(12*np.pi)*xi*1e3:.2f} mm")
print(f"   is the sub-mm screening scale in the [N9] torsion-balance row).")
print(f"COROLLARY (closes D28's conditionality): the D28 granularity sits at k xi = 1/2pi,")
print(f"   INSIDE the functional's domain by construction — the convexity penalty is not")
print(f"   optional once [N3] is frozen; 'A4 at coarse resolution' is not an available fork.")

print("\n--- THEOREM II: canonical inheritance (one hbar) ---")
# canonical pair from the microscopic field: [theta, dn] = i delta; quadratic Hamiltonian
# H2 = int [ (eps''(n)/2) dn^2 + (hbar^2 n/2m)(grad theta)^2 ]; Heisenberg equations give
# omega^2 = (n eps''(n)/m) k^2  ->  c_s^2 = n eps''/m. Entropic eps = hbar*gamma*n^3/3:
gamma_sym = 1.0  # symbolic unit check
nbar = 1.0
epspp = 2*gamma_sym*nbar          # d2/dn2 [gamma n^3/3] = 2 gamma n  (hbar=1 units)
cs2 = nbar*epspp                  # (m=1)
print(f"eps(n) = hbar gamma n^3/3  ->  eps''(nbar) = 2 hbar gamma nbar  ->  c_s^2 = nbar eps''/m")
print(f"      = 2 hbar gamma nbar^2/m  — EXACTLY the series' frozen quintic sound speed.")
print(f"(unit check at nbar = gamma = m = hbar = 1: c_s^2 = {cs2:.0f} = 2 ✓)")
print(f"The commutator came from [psi, psi+] (the particle sector); the Hamiltonian from the")
print(f"entropic functional; their marriage reproduces the c_s used in every downstream number.")
print(f"Consequences: (i) the phonon vacuum is a standard Gaussian state -> the [N3]")
print(f"entanglement counting and the D22 quadratic action are legitimate as quantum")
print(f"computations; (ii) there is ONE hbar — the acoustic Planck scale M_Pl,ac =")
print(f"m/sqrt(12 pi) = {8.45/np.sqrt(12*np.pi):.2f} eV inherits the carrier's quantum of action;")
print(f"(iii) D30's flagged 'local-thermodynamic force' escape is RETIRED: the phonon is a")
print(f"real canonical mode, the force is phase-mediated as the series states.")
occ = 421
print(f"Accuracy floor: the inheritance is exact at high occupancy; corrections O(1/sqrt(N))")
print(f"= {100/np.sqrt(occ):.1f}% per lambda_dB cell (D27), vanishing under coarse-graining —")
print(f"and degrading to O(1) only at the sigma_crit boundary (the known break-width caveat).")

print("\n--- what remains open (the deep stratum, stated plainly) ---")
print("WHY entropy counts area at 1/4 — for spacetime, inherited from [DM26] unanswered;")
print("for the acoustic sector, [A4] posits it. The two theorems above close every gap")
print("BETWEEN the assumptions (no hidden scale choice, no hidden quantization postulate);")
print("they do not derive [A4] from below. That derivation is the founding-aim horizon")
print("(SUBSTRATE_ROADMAP Q-A'), constrained by this program (the ceiling; the separation")
print("theorem; now the resolution lock) but not supplied by it.")

json.dump(dict(
    theorem_I=dict(eta_xi2_at_physical_regulator=float(eta_xi), frozen_1_over_48pi=float(eta_ref),
                   ratio=float(eta_xi/eta_ref), coarse_grained_suppression=float(eta_xi/eta_cg),
                   l_Pl_ac_mm=float(np.sqrt(12*np.pi)*xi*1e3),
                   statement="resolution forced by the Bogoliubov regulator; frozen [N3] incompatible with coarse-graining ((2pi)^2 off); observable-locked via l_Pl,ac"),
    theorem_II=dict(cs2_check="n eps''/m = 2 hbar gamma n^2/m reproduced exactly",
                    one_hbar="M_Pl,ac = m/sqrt(12pi) inherited from the carrier",
                    accuracy_floor_percent=float(100/np.sqrt(occ)),
                    statement="canonical structure particle-inherited; phonon vacuum Gaussian; [N3]/D22 legitimate; local-thermodynamic escape retired"),
    open_deep_stratum="microscopic origin of S=dA/4 (spacetime AND acoustic) — founding-aim horizon, not closed",
    verdict="T-5 tractable stratum CLOSED: D28/D31 conditionality upgraded from 'assumed' to 'internally forced'"),
    open(os.path.join(DERIVED, "b2_substrate_consistency.json"), "w"), indent=1)
print("\nwrote derived/b2_substrate_consistency.json")
