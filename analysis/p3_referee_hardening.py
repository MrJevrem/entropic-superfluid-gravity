#!/usr/bin/env python
"""
Paper III referee-hardening computations:
 (1) thermal-relic seduction: an 8.5 eV thermal boson gives the RIGHT abundance but
     a ~1e2 Mpc free-streaming length -> Lyman-alpha dead (sharpens 'thermal excluded')
 (2) occupancy Liouville pair: phase-space occupancy ~5e2 at cosmic mean AND in halos
     (coarse-grained phase-space density conserved) -> condensation criterion is
     virialization-proof
 (3) Mach numbers: Fornax GCs (Landau zero-drag claim) vs MW bar (v0.1 overclaim check)
 (4) wide-binary EFE regime: g_ext/a0 in the solar neighborhood
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import DERIVED

EV = 1.602176634e-19          # J
H = 6.62607015e-34
MSUN, PC = 1.98892e30, 3.0857e16
m_eV = 8.5
m_kg = m_eV * EV / (2.99792458e8) ** 2

print("=== (1) thermal seduction ===")
om_th = m_eV / 94.0     # neutrino-like decoupling
lam_fs = 31.0 * (30.0 / m_eV)   # Mpc, HDM free-streaming scaling
print(f"thermal 8.5 eV: Omega h^2 = {om_th:.3f} (measured 0.120 -- 'right-sized')")
print(f"but free-streaming ~ {lam_fs:.0f} Mpc -> hot dark matter, Lyman-alpha/structure dead")

print("\n=== (2) occupancy Liouville pair ===")
def occupancy(rho_kg_m3, v_ms):
    n = rho_kg_m3 / m_kg
    lam = H / (m_kg * v_ms)
    return n * lam ** 3
occ_cosmic = occupancy(0.26 * 9.47e-27, 3e3)          # mean DM density, linear-regime v~3 km/s
occ_halo = occupancy(0.01 * MSUN / PC ** 3, 2e5)      # 0.01 Msun/pc^3, v~200 km/s
print(f"cosmic mean (v=3 km/s):   n lambda^3 = {occ_cosmic:.0f}")
print(f"MW-like halo (v=200 km/s): n lambda^3 = {occ_halo:.0f}")
print(f"ratio {occ_halo/occ_cosmic:.2f} -- coarse-grained phase-space density ~conserved:")
print("  the condensation criterion survives virialization (Liouville), not by accident")

print("\n=== (3) Mach numbers (Landau-drag audit) ===")
for name, vf, vpert in (("Fornax GCs", 18.0, 10.0), ("MW bar @ corotation", 220.0, 220.0)):
    cs = vf / np.sqrt(2)
    M = vpert / cs
    verdict = "subsonic: ZERO drag (Landau)" if M < 1 else "SUPERSONIC: reduced-not-zero (phonon Cherenkov) -- fix scorecard row"
    print(f"{name:22s}: c_s = {cs:5.1f} km/s, v = {vpert:5.1f} -> Mach {M:.2f}  {verdict}")

print("\n=== (4) wide-binary EFE regime ===")
a0 = 1.080e-10
g_ext = 1.8e-10   # solar-neighborhood Galactic field
print(f"g_ext/a0 = {g_ext/a0:.2f} -> EFE-dominated regime; boost is EFE-suppressed by construction")
print("(the 0.12-0.17 band is a phonon-EFE computation, not AQUAL/QUMOND -- state formalism)")

json.dump(dict(omega_thermal_h2=om_th, lambda_fs_Mpc=lam_fs,
               occ_cosmic=occ_cosmic, occ_halo=occ_halo,
               mach_fornax=10 / (18 / np.sqrt(2)), mach_bar=np.sqrt(2),
               g_ext_over_a0=g_ext / a0),
          open(os.path.join(DERIVED, "p3_hardening.json"), "w"), indent=1)
print("\nsaved derived/p3_hardening.json")