#!/usr/bin/env python
"""
D30 — the condensation-rate computation (T-10 residue; user directive 2026-07-25).

Structure of the problem (from D29): condensation of the granulated halo field is
a TWO-STAGE process — (1) local kinetic quasi-condensation (speckle contrast decays
on microscopic times), then (2) PHASE ORDERING: the coherence length l_coh(t) must
grow to the scales the series' P(X) hydrodynamics actually uses (the phonon-EFE,
soliton, and R_SF computations all assume kpc-coherent theta). Stage 1 is fast
beyond doubt; the physical question is stage 2's growth law l_coh(t) ~ xi (t/tau)^beta.

This script: (a) the micro-timescales; (b) beta MEASURED from the D29 toy time
series; (c) the ladder of candidate ordering laws evaluated at the two physical
target scales (disc 30 kpc; R_SF 186 kpc); (d) the data-side requirement locks;
(e) the vortex inventory; (f) verdict.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import DERIVED

EV = 1.602176634e-19; C = 2.99792458e8; HBAR = 1.054571817e-34; H = 6.62607015e-34
PC = 3.0857e16; GYR = 3.156e16; AU = 1.496e11
m = 8.45*EV/C**2
vf = 156e3; sig = vf/np.sqrt(2); cs = vf/np.sqrt(2)
mu = 0.5*m*cs**2
tau_mic = HBAR/mu
xi = HBAR/(m*cs)
lam = H/(m*sig)
L_DISC = 30e3*PC; L_RSF = 186e3*PC
T_H = 13.8*GYR

print("=== D30: the condensation rate — from kinetics to phase ordering ===\n")

print("--- (1) micro-timescales (stage 1: local quasi-condensation) ---")
print(f"mu = m c_s^2/2 = {mu:.2e} J;  tau_micro = hbar/mu = {tau_mic*1e9:.2f} ns;  xi = {xi*1e3:.3f} mm")
print(f"local speckle-contrast decay: O(10-100) tau_micro ~ microseconds — INSTANT on every")
print(f"astrophysical clock. Stage 1 is never the bottleneck; the question is stage 2.")

print("\n--- (2) the ordering exponent measured in the D29 toy ---")
tb = json.load(open(os.path.join(DERIVED, "t10b_time_extension.json")))["series"]
t = np.array([p["t_tdyn"] for p in tb[3:]]); Cc = np.array([p["C"] for p in tb[3:]])
ell = np.sqrt(np.maximum(Cc, 1e-9))            # 2D: C ~ (l_coh/R)^2  ->  l ~ sqrt(C)
A = np.vstack([np.log(t), np.ones_like(t)]).T
beta, _ = np.linalg.lstsq(A, np.log(ell), rcond=None)[0]
print(f"t10b series (t/t_dyn = {t[0]:.0f}..{t[-1]:.0f}): l_coh ~ sqrt(C) fit  ->  beta = {beta:.2f}")
print(f"beta ~ 0.5-0.6 is the DIFFUSIVE coarsening class (2D vortex-pair annihilation,")
print(f"log-slow) — the known behavior of a T=0, dissipationless, 2D toy. Whether the")
print(f"physical system (3D line tangle + native normal-fraction dissipation + rotation)")
print(f"stays in this class or orders ballistically is exactly the fork below.")

print("\n--- (3) the ladder of ordering laws at physical scale ---")
print(f"{'law':34s}{'t(30 kpc)':>14}{'t(R_SF=186 kpc)':>18}")
rows = {}
for name, f in (
    ("ballistic (sound-limited), l/c_s", lambda L: L/cs),
    (f"toy exponent beta={beta:.2f}", lambda L: tau_mic*(L/xi)**(1/beta)),
    ("diffusive, tau*(l/xi)^2", lambda L: tau_mic*(L/xi)**2),
    ("Svistunov-type, l^2/(xi c_s)", lambda L: L**2/(xi*cs)),
):
    td, tr = f(L_DISC), f(L_RSF)
    fmt = lambda x: f"{x/GYR:10.2f} Gyr" if x < 1e3*GYR else f"10^{np.log10(x/GYR):.0f} Gyr"
    print(f"{name:34s}{fmt(td):>14}{fmt(tr):>18}")
    rows[name] = dict(t_disc_Gyr=float(td/GYR), t_RSF_Gyr=float(tr/GYR))
print("ballistic is the CAUSAL FLOOR (information moves at c_s); every diffusive-class law")
print(f"overshoots the age of the universe by 16-25 dex — there is no middle ground.")

print("\n--- (4) the data-side requirement locks ---")
t_ball_disc = L_DISC/cs; t_ball_rsf = L_RSF/cs
print(f"LOCK 1 (RAR tightness, ~0.13 dex): disc-scale coherence must be UNIFORM across")
print(f"  galaxies of different ages/histories -> tau(30 kpc) << spread of formation times ~ Gyr.")
print(f"  Ballistic gives {t_ball_disc/GYR:.2f} Gyr: PASS with margin. Diffusive-class: FAIL by >20 dex.")
print(f"LOCK 2 (a0(z) = const, favored 3-6 sigma): the force must be saturated by z ~ 1")
print(f"  -> tau(R_SF) < ~6 Gyr. Ballistic gives {t_ball_rsf/GYR:.2f} Gyr: PASS. Diffusive: FAIL.")
print(f"LOCK 3 (D22 consistency): the f_n(r) = (r/R_SF)^2 boundary zone WANTS lagging outer")
print(f"  coherence — ballistic ordering, which completes inside-out on ~Gyr, matches the shape.")

print("\n--- (5) the vortex inventory (what stage 2 must actually anneal) ---")
sp0 = lam; sp_eq = 2.9e9   # D23 equilibrium lattice spacing ~0.02 AU
print(f"granulated state: inter-vortex spacing ~ lambda_dB = {sp0*1e3:.1f} mm")
print(f"equilibrium (D23 rotating lattice): spacing ~ 0.02 AU = {sp_eq:.1e} m")
print(f"  -> ordering must coarsen the vortex spacing by {np.log10(sp_eq/sp0):.1f} dex —")
print(f"     the end state is NOT vortex-free (the halo rotates); it is the D23 lattice.")
print(f"The decisive physics: is 3D vortex-tangle decay in the quintic medium sound-limited")
print(f"(mutual friction on the native normal fraction; Kagan-Svistunov inside-out ordering)")
print(f"or diffusion-limited (T=0 pairwise annihilation, as the 2D toy shows)?")

print("\n--- (6) verdict ---")
print("The rate question is now an EXPONENT FORK, sharply posed:")
print("  SURVIVES  iff phase ordering is ballistic/sound-limited: tau = 0.27 Gyr (disc),")
print("            1.6 Gyr (R_SF) — passing all three data locks with margin.")
print("  DIES      if the toy's diffusive class persists to 3D physical scale: coherence")
print("            stalls ~16-25 dex short, the P(X) sector never turns on at kpc.")
print("The 2D T=0 toy sits in the diffusive class BY CONSTRUCTION (no dissipation channel,")
print("vortex-pair dynamics); the physical halo has the two ingredients known to accelerate")
print("ordering toward the sound-limited class: 3D line tension and mutual friction on the")
print("native normal component. NAMED DECISIVE COMPUTATION: 3D dissipative vortex-tangle")
print("decay in the quintic medium (folds into T-3, now with a pass/fail criterion).")
print("Flagged alternative (T-5-adjacent): if the entropic force is local-thermodynamic")
print("rather than phase-mediated, the kpc-coherence requirement itself relaxes — a")
print("formulation question upstream of this fork, not invoked here.")

json.dump(dict(mu_J=float(mu), tau_micro_ns=float(tau_mic*1e9), xi_mm=float(xi*1e3),
               beta_toy=float(beta), ladder=rows,
               locks=dict(rar_tightness="tau(30kpc) << Gyr", a0z_const="tau(R_SF) < ~6 Gyr",
                          d22_shape="inside-out ordering matches f_n(r)"),
               ballistic_Gyr=dict(disc=float(t_ball_disc/GYR), rsf=float(t_ball_rsf/GYR)),
               vortex_coarsening_dex=float(np.log10(sp_eq/sp0)),
               verdict="EXPONENT FORK: survives iff ordering is sound-limited (0.27/1.6 Gyr, all locks pass); "
                       "dies if the toy's diffusive class persists (16-25 dex short). Decisive: 3D dissipative "
                       "tangle decay in the quintic medium (-> T-3)."),
          open(os.path.join(DERIVED, "b2_wave_kinetic_rate.json"), "w"), indent=1)
print("\nwrote derived/b2_wave_kinetic_rate.json")
