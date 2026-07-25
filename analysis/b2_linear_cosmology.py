#!/usr/bin/env python
"""
D34 — linear cosmology beyond the four gates: branch selection + deviation ledger.

THE DISCOVERED ISSUE (stated first, resolved second): naively extrapolating the
in-halo quintic stiffness (fixed gamma, c_s ~ n) to the early universe gives
c_s(z_rec) > c — a relativistic-sound DM that erases structure and is excluded
instantly. Resolution — the stiffness is NOT a fixed microscopic gamma: D28
derived c_s = v_f/sqrt(2) as the EQUILIBRIUM property of the virialized
condensate (gamma_eff = m v^2/(4 hbar nbar^2), halo-specific). The entropic
response requires the horizon-fluctuation setting of a bound structure; the
homogeneous background has neither. Two self-consistent branches exist:
  DUST branch: background pressureless; only quantum pressure at the m-scale.
      Structure forms; halos virialize; the entropic sector activates inside
      them (D28) -> the observed universe. Selected by the misalignment initial
      condition (theta homogeneous and frozen: no gradients, no phonon sector).
  STIFF branch: entropic stiffness at cosmic mean -> w(z_rec) ~ O(0.1-0.3),
      no structure, no observers — and unreachable from the initial condition.
This is the same conclusion the frozen record reached as "shared-g Bogoliubov
corner FORBIDDEN (2.4 dex)"; here it is made a branch-selection statement.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import DERIVED

EV = 1.602176634e-19; C = 2.99792458e8; HBAR = 1.054571817e-34; G = 6.67430e-11
MPC = 3.0857e22; T_H = 4.35e17
m = 8.45*EV/C**2
rho_dm0 = 2.3e-27

print("=== D34: linear-cosmology branch selection and deviation ledger ===\n")
print("--- (1) the stiff branch, excluded ---")
cs_gal = 110e3; n_gal = 1e-22/m; n0 = rho_dm0/m
for z in (0, 1100):
    n = n0*(1+z)**3
    cs_naive = min(cs_gal*(n/n_gal), C/np.sqrt(2))
    w = (cs_naive/C)**2/3
    print(f"z = {z:5d}: fixed-gamma extrapolation c_s = {cs_naive:.2e} m/s -> w = {w:.2e}")
print("w(z_rec) ~ 0.17 is excluded by the CMB at ~inf sigma: the stiff branch is not the")
print("universe we are in — and it is UNREACHABLE: the misalignment field is born homogeneous")
print("and frozen (no gradients -> no phonon sector -> no stiffness to bootstrap). gamma is not")
print("fundamental: gamma_eff = m v_f^2/(4 hbar nbar^2) is the D28 equilibrium of a virialized")
print("halo. Background verdict: DUST (+ m-scale quantum pressure). [Matches the frozen gate.]\n")

print("--- (2) the deviation ledger on the dust branch ---")
lamJ = lambda rho: (HBAR**2/(G*rho*m**2))**0.25
for z, tag in ((0, "today"), (1100, "recombination")):
    rho = rho_dm0*(1+z)**3
    lj = lamJ(rho)
    kJ = 2*np.pi/lj*MPC   # in 1/Mpc comoving-ish (phys at z; conservative)
    print(f"{tag:14s}: quantum-pressure Jeans length = {lj:.2e} m = {lj/MPC:.1e} Mpc  ->  k_J ~ {kJ:.1e}/Mpc")
k_obs = 10.0
supp = (k_obs/ (2*np.pi/lamJ(rho_dm0)*MPC))**4
print(f"P(k) suppression at k = {k_obs}/Mpc: (k/k_J)^4 ~ {supp:.1e} — {abs(np.log10(supp)):.0f} dex below")
print(f"observability; CMB and linear P(k) are CDM to every measurable digit.\n")

print("--- (3) the a0-channel at recombination: two gates ---")
# gate A: the D33 tail — typical peak-scale accelerations vs a0
A0 = 1.080e-10
for k in (0.01, 0.1, 1.0):
    lam_phys = (1/k)*MPC/1101
    rho_m = 9.9e-27*1101**3
    g_pert = 4*np.pi*G*rho_m*1e-5*lam_phys/3
    x = g_pert/A0
    print(f"k = {k:5.2f}/Mpc: g_pert(z=1100) ~ {g_pert:.1e} m/s^2, x = g/a0 = {x:8.2f} -> tail factor e^(-x/pi) = {np.exp(-min(x,500)/np.pi):.2f}")
print("gate A alone is NOT sufficient at all scales (x ~ O(1-100) across the peaks) — the")
print("operative gate is B: ESTABLISHMENT. The rectified channel is a property of the")
print("condensate's acoustic-horizon setting; on the dust branch the medium's collective")
print("speed at these scales is the quantum speed c_q = hbar k/(2m):")
for k in (0.01, 0.1, 1.0):
    lam_phys = (1/k)*MPC/1101
    cq = HBAR*(2*np.pi/lam_phys)/(2*m)
    t_est = lam_phys/cq
    print(f"k = {k:5.2f}/Mpc: c_q = {cq:.1e} m/s; establishment time lambda/c_q = {t_est/T_H:.1e} t_H")
print("=> the channel cannot form across linear-cosmology scales by 20+ dex in time — the")
print("frozen 'sound-crossing protection' gate, recomputed and extended to a full k-ladder.\n")

print("--- (4) verdict and the named residue ---")
print("Linear cosmology = CDM is now a THEOREM-SHAPED statement on the dust branch: quantum")
print("pressure 11-12 dex beyond reach; the entropic channel doubly gated (tail + establishment);")
print("the stiff branch excluded and unreachable. Residues (production): a full Boltzmann run")
print("(cosmetic given the ledger) and the nonlinear/halo-scale wave simulations (real).")

json.dump(dict(stiff_branch=dict(w_rec=0.17, status="excluded + unreachable (misalignment IC)"),
               dust_ledger=dict(kJ_today_per_Mpc=float(2*np.pi/lamJ(rho_dm0)*MPC),
                                Pk_suppression_at_k10=float(supp)),
               a0_channel=dict(gateA="x ~ O(1-100) at peaks — insufficient alone",
                               gateB="establishment time 1e20+ t_H — operative"),
               verdict="linear cosmology = CDM as branch-selection theorem + deviation ledger; residue: Boltzmann run (cosmetic), nonlinear wave sims (real)"),
          open(os.path.join(DERIVED, "b2_linear_cosmology.json"), "w"), indent=1)
print("\nwrote derived/b2_linear_cosmology.json")
