#!/usr/bin/env python
"""
D36 / T-2 — the supersonic drag computation (the bar), done honestly.

The [N9] scorecard row asserted "reduced-not-zero drag on the Mach-1.4 bar
(open computation, flagged)". This computes it. Linear response of a T=0
superfluid to a supersonic gravitational perturber = the compressible-medium
(Ostriker) result: I_sf(M) = (1/2) ln(1 - 1/M^2) + ln(Lambda) for M > 1,
against the collisionless (Chandrasekhar) I_cdm = ln(Lambda) [erf(X) -
2X e^{-X^2}/sqrt(pi)], X = v/(sqrt(2) sigma). Subsonic superfluid: zero
steady-state drag (Landau) — the Fornax side. Quintic dispersion correction:
Cherenkov band limited to k xi < sqrt(2(M^2-1)) — at bar scales k xi ~ 1e-24,
irrelevant. Residual-tangle mutual friction (D31): f_n-suppressed, negligible.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.special import erf
from t3_guard import DERIVED

G = 6.67430e-11; MSUN = 1.98892e30; PC = 3.0857e16; GYR = 3.156e16
cs = 110e3; sigma = 110e3
v_bar = 156e3                      # bar-halo relative velocity scale (pattern-speed regime)
Mmach = v_bar/cs
lnL = np.log(3e3/0.3e3)*0 + np.log(10.0)   # r_max/r_min ~ 3 kpc / 0.3 kpc -> ln 10
M_bar = 1e10*MSUN; rho = 0.01*MSUN/PC**3

print("=== D36: supersonic drag on the Galactic bar — the honest numbers ===\n")
I_sf = 0.5*np.log(1 - 1/Mmach**2) + lnL
X = v_bar/(np.sqrt(2)*sigma)
I_cdm = lnL*(erf(X) - 2*X*np.exp(-X**2)/np.sqrt(np.pi))
print(f"Mach = {Mmach:.2f} (v = {v_bar/1e3:.0f} km/s vs c_s = {cs/1e3:.0f});  ln Lambda = {lnL:.2f}")
print(f"I_superfluid (Ostriker, M>1) = {I_sf:.2f}")
print(f"I_CDM (Chandrasekhar, X={X:.2f}) = {I_cdm:.2f}")
print(f"ratio I_sf/I_cdm = {I_sf/I_cdm:.2f}")
F = 4*np.pi*G**2*M_bar**2*rho/v_bar**2
tau_sf = M_bar*v_bar/( F*I_sf )
tau_cdm = M_bar*v_bar/( F*I_cdm )
print(f"point-mass slowdown timescale (ORDER-OF-MAGNITUDE ONLY — real bar torque is resonant):")
print(f"  superfluid {tau_sf/GYR:.2f} Gyr vs CDM {tau_cdm/GYR:.2f} Gyr — the ROBUST deliverable is the ratio.\n")

print("--- verdict, stated adversely where it is adverse ---")
print(f"AT FACE VALUE the phonon-radiation drag at Mach {Mmach:.1f} is TWICE the collisionless")
print(f"drag (ratio {I_sf/I_cdm:.2f}) — near Mach 1 the compressible wake is MORE efficient than the")
print(f"collisionless one. The earlier scorecard wording 'eased' is NOT supported: it is inverted.")
print(f"The fast-bar tension afflicts this framework at least as hard as CDM. RECORDED AS ADVERSE. Mitigations that are")
print(f"real but unquantified here: (i) the bar is an extended rotating quadrupole, not a point")
print(f"mass — the resonant (Weinberg-class) torque differs from Ostriker by O(1) factors either")
print(f"way; (ii) partial corotation of the inner condensate (vortex-lattice spin-up) reduces")
print(f"the effective Mach number toward the zero-drag boundary; (iii) near-sonic linear theory")
print(f"is marginal at M = 1.4. The honest status: OPEN TENSION, shared with CDM, not resolved")
print(f"by superfluidity at this Mach number; the [N9] bar row must be reworded at the next")
print(f"paper pass.\n")

print("--- the subsonic side (Fornax) — robust ---")
Mach_f = 0.8
print(f"Fornax GCs at Mach {Mach_f}: subsonic -> ZERO steady-state linear drag (Landau criterion;")
print(f"the transient wake decays). The classic Fornax timing problem is absent by structure —")
print(f"this half of the scorecard row stands, and is the discriminant against CDM [M].")
print(f"Quintic Cherenkov band: k xi < sqrt(2(M^2-1)) = {np.sqrt(2*(Mmach**2-1)):.2f} — at bar scales")
print(f"k xi ~ 1e-24: no additional gating. Residual-tangle mutual friction: f_n-suppressed, ~1e-2")
print(f"of the phonon drag at most: negligible either way.")

json.dump(dict(mach=float(Mmach), I_sf=float(I_sf), I_cdm=float(I_cdm), ratio=float(I_sf/I_cdm),
               tau_sf_Gyr=float(tau_sf/GYR), tau_cdm_Gyr=float(tau_cdm/GYR),
               verdict="ADVERSE: Mach-1.4 phonon-radiation drag is ~2x the collisionless drag; 'eased' wording "
                       "unsupported; fast-bar tension shared with CDM, unresolved; mitigations named "
                       "(resonant torque; partial corotation). Fornax subsonic zero-drag ROBUST [M].",
               paper_action="reword [N9] bar row at next pass"),
          open(os.path.join(DERIVED, "b2_supersonic_drag.json"), "w"), indent=1)
print("\nwrote derived/b2_supersonic_drag.json")
