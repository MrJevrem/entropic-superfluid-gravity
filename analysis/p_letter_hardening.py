#!/usr/bin/env python
"""
Letter referee-hardening: verify every quantitative claim, with error propagation.
 (1) a0_pred = c^2 sqrt(Lambda/24pi) with Planck error bar (currently absent from draft)
 (2) tension vs BOTH comparators: canonical RAR g+ = 1.20 +/- 0.02 +/- 0.24 (McGaugh+16 PRL)
     and marginalized 1.128 +/- 0.019 (nuisance-marginalized fit) -- draft quotes only 1.128
 (3) modular candidate c^2 sqrt(L/3)/2pi and the exact sqrt(pi/2) ratio
 (4) epsilon ledger at PLANCK Omega_L=0.6847 (draft mixed 0.685 and 0.70)
 (5) implied H0: at Omega_L=0.6847 vs 0.70 -- hygiene check on the '69.6' corollary
 (6) Lambda-meter: % vs Planck for both comparators + systematic on Lambda
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import DERIVED

C = 2.99792458e8
MPC = 3.0856775814913673e22
OL, dOL = 0.6847, 0.0073          # Planck 2018 TT,TE,EE+lowE+lensing
H0P, dH0P = 67.36, 0.54           # km/s/Mpc
h0 = H0P * 1e3 / MPC
LAM = 3 * OL * h0**2 / C**2
dLAM = LAM * np.sqrt((dOL/OL)**2 + (2*dH0P/H0P)**2)

a0_pred = C**2 * np.sqrt(LAM / (24*np.pi))
da0_pred = a0_pred * 0.5 * dLAM/LAM
a0_mod = C**2 * np.sqrt(LAM/3) / (2*np.pi)

comparators = {"canonical RAR (McGaugh+16)": (1.20e-10, 0.02e-10, 0.24e-10),
               "marginalized fit":           (1.128e-10, 0.019e-10, 0.226e-10)}  # 20% M/L syst

print("=== Letter hardening ===")
print(f"Lambda(Planck) = {LAM:.4e} m^-2  (+/- {100*dLAM/LAM:.1f}%)")
print(f"a0_pred  = ({a0_pred*1e10:.3f} +/- {da0_pred*1e10:.3f}) e-10   <- prediction now HAS an error bar")
print(f"a0_mod   = {a0_mod*1e10:.3f} e-10 ; pred/mod = {a0_pred/a0_mod:.4f}  vs sqrt(pi/2) = {np.sqrt(np.pi/2):.4f}")
res = {"a0_pred": a0_pred, "da0_pred": da0_pred, "a0_modular": a0_mod, "Lambda": LAM}
for name, (v, ds, dsys) in comparators.items():
    gap = (v - a0_pred)/a0_pred
    sig = (v - a0_pred)/np.sqrt(ds**2 + dsys**2 + da0_pred**2)
    print(f"  vs {name:28s}: obs {v*1e10:.3f}  pred low by {100*gap:.1f}%  ({sig:+.2f} sigma incl. syst)")
    res[name] = dict(obs=v, gap_pct=100*gap, sigma=float(sig))
    lam_meter = 24*np.pi*v**2/C**4
    print(f"     Lambda-meter: {lam_meter:.3e}  = Planck {100*(lam_meter/LAM-1):+.0f}%  (syst on Lambda ~{200*np.sqrt(ds**2+dsys**2)/v:.0f}%)")

print("\n--- epsilon ledger (ALL at Planck Omega_L = 0.6847) ---")
eps_pred = np.sqrt(OL/(8*np.pi))
print(f"eps_pred = sqrt(OL/8pi) = {eps_pred:.4f} ; 1/6 = {1/6:.4f} (diff {100*(1/6/eps_pred-1):+.2f}%) ; "
      f"1/2pi = {1/(2*np.pi):.4f} (diff {100*(1/(2*np.pi)/eps_pred-1):+.2f}%)")
for H in (67.36, 69.6, 70.3, 73.04):
    ch = C * H*1e3/MPC
    print(f"  eps_meas(a0=1.128, H0={H}) = {1.128e-10/ch:.4f}")

print("\n--- implied H0 hygiene check (a0 -> Lambda -> H0 at fixed Omega_L) ---")
for a0v, tag in ((1.128e-10, "marginalized"), (1.20e-10, "canonical")):
    lam_a0 = 24*np.pi*a0v**2/C**4
    for ol in (0.6847, 0.70):
        H0i = np.sqrt(lam_a0*C**2/(3*ol)) * MPC/1e3
        print(f"  a0={a0v*1e10:.3f} ({tag:12s}), Omega_L={ol}:  H0 = {H0i:.1f}")
res["H0_implied_planckOL"] = float(np.sqrt(24*np.pi*(1.128e-10)**2/C**2/(3*OL))*MPC/1e3)
json.dump({k: (v if not isinstance(v, dict) else v) for k, v in res.items()},
          open(os.path.join(DERIVED, "p_letter_hardening.json"), "w"), indent=1, default=float)
print("\nsaved derived/p_letter_hardening.json")