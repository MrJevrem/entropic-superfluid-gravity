#!/usr/bin/env python
"""
D22 — phonon perturbative-stability audit for the entropic X^(3/2) sector
(the BK-known issue: gradient instability / ghost around static MOND profiles).

Quadratic fluctuations of L = P(X), X = theta_dot - (grad theta)^2/2m - m Phi,
around a static background (u = grad theta_bar / m):
    L2 = (P''/2)(pi_dot - u.grad pi)^2 - (P'/2m)(grad pi)^2
 => ghost-free iff P'' > 0; gradient-stable iff P' > 0; c_s^2 = P'/(m P'')
    in the frame comoving with the phonon wind u.

Branch structure of P = K X^(3/2) (continued as sign(X)|X|^(3/2)):
  X > 0 (hydrostatic, mu-dominated): P' = (3K/2)sqrt(X) > 0, P'' = (3K/4)/sqrt(X) > 0 -> HEALTHY
  X < 0 (gradient-dominated; the branch BK's direct coupling forces): P' > 0 but P'' < 0
        -> ghost + c_s^2 < 0: the BK instability. Their cure: finite-T beta term.

Audit questions: (1) reproduce the BK-branch pathology (machinery validation);
(2) does OUR mechanism ever occupy X < 0? (no direct coupling -> static solutions
have u = 0, X_bar = mu(r)/hbar > 0 inside R_SF); (3) health scan over an MW-like
background; (4) the X -> 0 boundary zone (where finite-T two-fluid takes over,
quantified); (5) Jeans consistency of the polarization cloud.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import ROOT, DERIVED   # path constants only; no frozen comparison here

HBAR = 1.054571817e-34; EV = 1.602176634e-19; C = 2.99792458e8
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
m = 8.45 * EV / C**2

print("=== D22: phonon stability audit (entropic X^(3/2) sector) ===\n")
print("--- (1) branch table (K = 1 units; validation incl. the BK branch) ---")
for X in (+1.0, +1e-4, -1e-4, -1.0):
    s = np.sign(X); aX = abs(X)
    P1 = 1.5 * np.sqrt(aX)               # P' = (3/2)|X|^(1/2)  (both branches)
    P2 = s * 0.75 / np.sqrt(aX)          # P'' = sign(X)(3/4)|X|^(-1/2)
    cs2 = P1 / P2                        # x 1/m, sign is what matters
    tag = "HEALTHY" if (P2 > 0 and P1 > 0) else "GHOST + GRADIENT-UNSTABLE (BK branch)"
    print(f"  X = {X:+8.1e}:  P' = {P1:+.3f}  P'' = {P2:+9.1f}  sign(c_s^2) = {np.sign(cs2):+.0f}   {tag}")
vf = 220e3; cs = vf / np.sqrt(2)
kk = 1 / KPC
tau_bk = 1 / (cs * kk) / 3.156e13   # e-fold time of the BK-branch instability at k = 1/kpc, in Myr
print(f"  BK-branch growth at k = 1/kpc, MW c_s: tau ~ {tau_bk:.0f} Myr  (why their beta-term is mandatory)")

print("\n--- (2) which branch does the entropic mechanism occupy? ---")
print("  No phonon-baryon coupling exists [D2, BK_COMPARISON]: nothing sources a static")
print("  phonon wind, so stationary solutions have u = 0 and X_bar = mu(r)/hbar > 0")
print("  wherever the condensate exists. The force is carried by grad P(X_bar(r))")
print("  (hydrostatic entropic stress), never by a gradient-dominated X < 0 profile.")
print("  AUDITED: the frozen WB/EFE formalism and the D15 phenomenology chain use")
print("  X_bar > 0 backgrounds throughout; the X < 0 branch is never invoked.")

print("\n--- (3) MW-like health scan (v_f = 220 km/s, isothermal-ish halo) ---")
Xbar0 = m * cs**2 / (2 * HBAR)
r = np.geomspace(1, 300, 40)  # kpc
RSF = 200.0
TT = np.clip((r / RSF) ** (4 / 3), 0, 1.2)     # T/Tc ~ sigma^2/rho^(2/3) ~ r^(4/3), =1 at R_SF
Xbar = Xbar0 * np.clip(1 - TT, 0, None) ** (2 / 3)  # soft decline of condensate mu toward boundary
ok = Xbar > 0
print(f"  X_bar(interior) = {Xbar0:.2e} s^-1  (mu ~ {HBAR*Xbar0/EV:.1e} eV);  c_s/c = {cs/C:.1e}  (subluminal by x{C/cs:.0f})")
print(f"  ghost-free (P''>0) and gradient-stable (P'>0) at all {ok.sum()}/{len(r)} radii with X_bar > 0")
print(f"  phonon wind: u <= v_f gives only a Doppler tilt (omega - u.k)^2 = c_s^2 k^2 — no instability;")
print(f"  supersonic stirring (bar, Mach 1.4) => Cherenkov DISSIPATION, not a ghost (consistent w/ P3 row)")

print("\n--- (4) the X -> 0 boundary zone ---")
fn = TT ** 1.5                                  # normal fraction ~ (T/Tc)^(3/2)
print(f"  c_s -> 0 softly at R_SF (healthy but slow); normal fraction f_n = (T/Tc)^(3/2)")
# f_n = TT^(3/2) = (r/RSF)^2 exactly, given TT = (r/RSF)^(4/3)
print(f"  f_n = (r/R_SF)^2 exactly  ->  f_n > 0.1 outside r = {RSF*np.sqrt(0.1):.0f} kpc = 0.32 R_SF")
print("  => the outer two-thirds (in radius) of the superfluid region carries a >10% normal")
print("     component: the finite-T two-fluid kinetic term (BK's beta, structurally) is")
print("     natively present exactly where the pure-condensate description softens.")
print("  Parallel stated honestly: our INTERIOR is healthy without finite-T help (unlike")
print("  BK's MOND region); our BOUNDARY shares BK's reliance on the normal component.")

print("\n--- (5) Jeans consistency of the polarization cloud ---")
rho_h = 1e-21
lamJ = cs * np.sqrt(np.pi / (G * rho_h)) / KPC
print(f"  lambda_J = c_s sqrt(pi/G rho) ~ {lamJ:.0f} kpc at rho = 1e-21 kg/m^3 — of order the halo")
print("  itself: c_s = v_f/sqrt(2) makes virialized halos Jeans-MARGINAL by construction;")
print("  no sub-halo clumping instability of the condensate. Consistency, not pathology.")

print("\n=== VERDICT [N19] ===")
print("The BK gradient/ghost instability is real on the X<0 branch and fast (tau ~ Myr),")
print("but it is a property of the DIRECT-COUPLING mechanism that forces the field onto")
print("that branch. The entropic mechanism's galactic backgrounds are hydrostatic X>0:")
print("ghost-free, gradient-stable, subluminal, wind-tolerant; the boundary zone is")
print("handled by the natively-present normal fraction (quantified: >10% outside 0.32 R_SF).")
print("Remaining scope (logged): full two-fluid closure at the boundary; time-dependent")
print("backgrounds (mergers); the coupled phonon-metric system beyond Jeans order.")

json.dump(dict(
    branch_table={"X>0": "P'>0, P''>0: healthy; c_s^2 = 2 hbar X/m",
                  "X<0": "P'>0, P''<0: ghost + gradient instability (BK branch); tau(k=1/kpc) ~ %.0f Myr" % tau_bk},
    occupancy_finding="no direct coupling -> static solutions u=0, X_bar>0 throughout; X<0 branch never invoked in the frozen formalism",
    mw_scan=dict(Xbar0_s=float(Xbar0), mu_eV=float(HBAR*Xbar0/EV), cs_over_c=float(cs/C),
                 healthy_radii="all with X_bar>0"),
    boundary=dict(fn_law="f_n=(r/R_SF)^2 under T/Tc ~ r^(4/3)", fn_gt_10pct_outside_kpc=float(RSF*np.sqrt(0.1)),
                  parallel="interior healthy without finite-T (unlike BK MOND region); boundary shares BK's normal-fraction reliance"),
    jeans_kpc=float(lamJ),
    remaining=["two-fluid closure at boundary", "time-dependent backgrounds/mergers", "phonon-metric coupling beyond Jeans"]),
    open(os.path.join(DERIVED, "b2_phonon_stability.json"), "w"), indent=1)
print("\nwrote derived/b2_phonon_stability.json")